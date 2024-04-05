// mpirun -np 64 -bind-to core:2 ./surf_pdhg -N 4 -m data/dino0.vtk -cgI 100 -tC 4 -alpha 0.1 -s -or 0 
// mpirun -np 64 -bind-to core:2 ./surf_pdhg -N 3 -m data/shell_quad.mesh -cgI 50 -tC 2 -alpha 0.0 -or 2

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "kron_quad.hpp"
#include "wass_laplace.hpp"
#include "wass_multigrid.hpp"
#include "bary_rhs.hpp"

// optimization solver
#include "brent.hpp"

using namespace std;
using namespace mfem;

double tolPDHG = 1e-3; // PDHG tolerance
double sigma_u = 1.0, sigma_phi = 1.0; // PDHG parameter: global parameter
int iterPnt = 100; // print every iterPnt steps

int N = 2; // N : # of species, N: # of reaction rates
// tunable GLOBAL parameters
// Diffusion parameters
double beta = 0.0;// interaction rate (set to zero by default)

// V1: always taken to be V1(rho) = rho
// V2 = alpha *(rho1-rho2)/(log(rho1)-log(rho2))
double alpha = 0.0; 
int testCase = 1; // test case 1: gaussian translation
                  // test case 2: image transfer
int cgI=10, ptI = 100; // maximum cg iteration & max brent iteration

// brent solver parameters
double rMin = 1e-6, rMax = 40.0; 
const int double_bits = std::numeric_limits<double>::digits / 2;

double rhoInit(const Vector &, double, int );

double   V1(double);
double   V2(double, double, double);
double   E(double);

// optimization function
double  F(double rho, Vector rhos,  double rhobar, 
    double mbar2, Vector sbars, int k);

// compute L1norm of a quadrature function
double L1Norm(QuadratureFunction qf);

// scale mesh
bool scale=false;

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();
   
   const char *mesh_file = "data/shell_quad.mesh";
   
   int n_time = 4;
   int ser_ref_levels = 0;
   int par_ref_levels = 0;
   int h_refs = 0;
   int p_refs = 0;
   bool paraview = true;
   int iterALGMax = 10000;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&n_time, "-nt", "--n-coarse-time",
                  "Number of coarse time elements.");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&iterALGMax, "-alg", "--iterALG",
                  "Total ALG iteration counts");
   args.AddOption(&tolPDHG, "-tolPDHG", "--tolPDHG",
                  "PDHG Tolerance");
   args.AddOption(&h_refs, "-gr", "--geometric-refinements",
                  "Number of geometric refinements done prior to order refinements.");
   args.AddOption(&p_refs, "-or", "--order-refinements",
                  "Number of order refinements. Finest level in the hierarchy has order 2^{or}.");
   args.AddOption(&paraview, "-pv", "--paraview", "-no-pv", "--no-paraview",
                  "Save data files for ParaView visualization.");
   args.AddOption(&iterPnt, "-iP", "--iterPnt", "print every # steps");
   args.AddOption(&scale, "-s", "--scale", "-no-s", "--no-scale",
                  "Scaling mesh");
   // interaction & reaction
   args.AddOption(&N, "-N", "--N", "# of species");
   args.AddOption(&alpha, "-alpha", "--alpha", "reaction strength");
   args.AddOption(&beta, "-beta", "--beta", "interaction rate");
   args.AddOption(&testCase, "-tC", "--tC", "testCase");
   args.AddOption(&rMax, "-rMax", "--rMax", "max density val in brent");
   args.AddOption(&cgI, "-cgI", "--cgI", "# cg iteration");
   args.AddOption(&ptI, "-ptI", "--ptI", "# brent iteration");
   args.ParseCheck();
   
   ParMesh coarse_space_mesh = [&]()
   {
      // READ from an mfem-compatible mesh
      Mesh serial_mesh(mesh_file);
      if (scale){
         double factor = 0.01;
         GridFunction *nodes = serial_mesh.GetNodes();
         if (nodes == NULL)
         {
            for (int i = 0; i < serial_mesh.GetNV(); i++)
            {
               double *v = serial_mesh.GetVertex(i);
               for (int d = 0; d < serial_mesh.SpaceDimension(); d++)
                 v[d] *= factor;
            }
         }
         else
         {
            *nodes *= factor;
         }
      }

      for (int l = 0; l < ser_ref_levels; l++) { serial_mesh.UniformRefinement(); }
      ParMesh par_mesh(MPI_COMM_WORLD, serial_mesh);
      serial_mesh.Clear();
      for (int l = 0; l < par_ref_levels; l++) { par_mesh.UniformRefinement(); }
      return par_mesh;
   }();

   ParMesh coarse_time_mesh = [&]()
   {
      MPI_Group world_group;
      MPI_Comm_group(MPI_COMM_WORLD, &world_group);
      MPI_Group local_group;
      int ranks[1] = { Mpi::WorldRank() };
      MPI_Group_incl(world_group, 1, ranks, &local_group);
      MPI_Comm local_comm;
      MPI_Comm_create(MPI_COMM_WORLD, local_group, &local_comm);
      Mesh serial_mesh = Mesh::MakeCartesian1D(n_time);
      return ParMesh(local_comm, serial_mesh);
   }();

   MultigridHierarchy space_hierarchy(coarse_space_mesh, h_refs, p_refs);
   MultigridHierarchy time_hierarchy(coarse_time_mesh, h_refs, p_refs);

   const int order_reduction = 1;
   const int order_h1 = space_hierarchy.GetOrder();
   const int order_l2 = order_h1 - order_reduction;
   ParMesh &space_mesh = space_hierarchy.GetFinestMesh();
   ParMesh &time_mesh = time_hierarchy.GetFinestMesh();
   ParFiniteElementSpace &Vs = space_hierarchy.GetFinestSpace();
   ParFiniteElementSpace &Vt = time_hierarchy.GetFinestSpace();

   const int dim_s = space_mesh.SpaceDimension(); // background space dimension
   const int dim_sf = space_mesh.Dimension(); // surface  dimension
   const int dim = dim_s + 1; // spacetime dimension
   
   const int nVs = Vs.GetTrueVSize();
   const int nVt = Vt.GetTrueVSize();
   const HYPRE_BigInt nVs_global = Vs.GlobalTrueVSize();
    
    string folder = "ParaView/surface/";
   
    
    folder += "TestCase"+to_string(testCase);
    folder += "/" + to_string(Vs.GlobalTrueVSize()) + 'x' + to_string(nVt) + '/';
    if (Mpi::Root()){
      string cmd = "mkdir -p " + folder;
      int mkdir = system(cmd.c_str());
    }
    string file = folder + "err.txt";
   {
      
      if (Mpi::Root())
      {
        cout << "Number of space unknowns: " << nVs_global << '\n'
             << "Number of time unknowns:  " << nVt << '\n'
             << "Total number of unknowns: " << nVs_global*nVt << endl;
      }
   }

   GridFunction time_nodes(&Vt); // time coordinates
   time_mesh.GetNodes(time_nodes);

   const Geometry::Type geom_s = space_mesh.GetElementGeometry(0);
   const Geometry::Type geom_t = time_mesh.GetElementGeometry(0);
   const int q_gl = 2*order_l2 + 1;

   const IntegrationRule &ir_s = IntRules.Get(geom_s, q_gl);
   const IntegrationRule &ir_t = IntRules.Get(geom_t, q_gl);

   QuadratureSpace Qs(space_mesh, ir_s);
   QuadratureSpace Qt(time_mesh, ir_t);

   const double mass_coeff = 2.0;
   KroneckerLaplacian kron(Vs, Vt, mass_coeff, n_time);
   const int nV = kron.Height();

   KroneckerMultigrid mg(space_hierarchy.GetSpaceHierarchy(),
                         time_hierarchy.GetSpaceHierarchy(),
                         mass_coeff,
                         n_time);

   CGSolver cg(MPI_COMM_WORLD);
   cg.SetPrintLevel(IterativeSolver::PrintLevel().None());
   cg.SetOperator(kron);
   cg.SetPreconditioner(mg);
   cg.SetMaxIter(cgI);
   cg.SetRelTol(1e-4);
   cg.SetAbsTol(1e-6);

   DenseMatrix X(nV, N), dX(nV, N);
   Vector X0(nV), dX0(nV);
   X = 0.0; dX = 0.0; 

   // (m & rho)
   KroneckerQuadratureFunction u_qf(Qs,  Qt, dim*N); // (mx, my, rho)
   //KroneckerQuadratureFunction u_qf1(Qs, Qt, dim*N); // old data
   KroneckerQuadratureFunction du_qf(Qs, Qt, dim*N); // grad(phi)
   
   // source terms
   KroneckerQuadratureFunction s_qf(Qs,  Qt, N); // source
   //KroneckerQuadratureFunction s_qf1(Qs, Qt, N); // old data
   KroneckerQuadratureFunction ds_qf(Qs, Qt, N);  // phi_gf
   
   // phi @ quad pts
   KroneckerQuadratureFunction phi_qf(Qs, Qt, N);      // phi_gf
   KroneckerQuadratureFunction dphi_qf0(Qs, Qt, dim);  // grad(phi_gf)
   KroneckerQuadratureFunction phi_qf0(Qs, Qt, 1);     // phi_gf
   
   // terminal data 
   KroneckerFaceInterpolator face_interp(Vs, Qs);
   QuadratureFunction rhoT_qf(Qs);
   QuadratureFunction rhoT_qf1(Qs);
   QuadratureFunction drhoT_qf(Qs);
   QuadratureFunction phiT_qf(Qs,N);
   QuadratureFunction phiT_qf0(Qs);

   // For RHS vector in Phi Solver
   QuadratureFunction rhsT_qf(Qs); 
   QuadratureFunctionCoefficient rhsT_coeff(rhsT_qf);
   
   const int nW = Qs.GetSize()*Qt.GetSize();
   const int nWf = Qs.GetSize();

   u_qf = 0.0;
   du_qf = 0.0;
   s_qf = 0.0;
   ds_qf = 0.0;
   rhoT_qf = 0.0;
   phi_qf = 0.0;
   phiT_qf = 0.0;
   
   // incremental
   //u_qf1 = 0.0;
   //s_qf1 = 0.0;
   rhoT_qf1 = 0.0;

   // INITIAL GUESS for VOL (u_qf)  and BND (rhoT_qf)
   for (int k=0; k<N; k++)
   {
      KroneckerQuadratureFunction rho(Qs, Qt, 1);
      auto rho0 = [k](const Vector &x, double t) {return rhoInit(x, t, k);};
      FunctionCoefficient rho_coeff(rho0);
      rho.Project(rho_coeff);
      for (int i = 0; i < rho.Size(); ++i)
         u_qf[i*dim*N+k*dim+dim-1] = rho[i];
      
      QuadratureFunction rhoT(Qs);
      auto rhoT0 = [k](const Vector &x) {return rhoInit(x, 0, k);};
      FunctionCoefficient rhoT_coeff(rhoT0);
      rhoT_coeff.Project(rhoT);
      for (int i = 0; i < rhoT.Size(); ++i)
         rhoT_qf[i] += rhoT[i];
   }
   // BND/VOL term get average
   rhoT_qf /= double(N);

   // Right-hand side
   KroneckerLinearForm B(Vs, Vt, Qs, Qt);

   DenseMatrix B_bdr(nVs, N);
   Vector B_bdr_0(nVs),  B_bdr_T(nVs);
   
   // INITIAL DATA ON BOTTOM
   // Bottom boundary (t = 0)
   for (int k=0;k<N;k++){
     auto rho0 = [k](const Vector &x) {return -rhoInit(x, 0.0, k);};
     FunctionCoefficient rho0_coeff(rho0);
     
     ParLinearForm b_bdr_0(&Vs);
     b_bdr_0.AddDomainIntegrator(new DomainLFIntegrator(rho0_coeff));
     for (auto *i : *b_bdr_0.GetDLFI()) { i->SetIntRule(&ir_s); }
     b_bdr_0.Assemble();
     b_bdr_0.ParallelAssemble(B_bdr_0);
     B_bdr.SetCol(k, B_bdr_0);
   }

   ParLinearForm b_bdr_T(&Vs);
   b_bdr_T.AddDomainIntegrator(new DomainLFIntegrator(rhsT_coeff));
   for (auto *i : *b_bdr_T.GetDLFI()) { i->SetIntRule(&ir_s); }
   b_bdr_T.UseFastAssembly(true);
   
   // top bdry updates
   auto set_rhs_qf = [&](int k)
   {
      // the face dofs (on top bdry)
      for (int i=0; i< nWf; i++){
         rhsT_qf[i] = sigma_phi * rhoT_qf[i];
         for (int k2=0; k2<k; k2++){
            rhsT_qf[i] -= phiT_qf[i*N+k2];
         }
      }
   };

   if (Mpi::Root())
   {
      std::cout << std::string(75, '=') << std::endl;
   }

   StopWatch sw_a, sw_b;
   double rt[2];
   
   int cgIt[N];
   double cgErr[N];
   
   // visualization
   auto vis = [&](bool lor)
   {
      ParMesh lor_space_mesh;
      if (lor)
      {
         lor_space_mesh = ParMesh::MakeRefined(space_mesh, order_h1, Quadrature1D::GaussLobatto);
      }
      ParMesh &vis_mesh = lor ? lor_space_mesh : space_mesh;
      const int vis_order = lor ? 0 : order_l2;

      L2_FECollection l2_space_fec(vis_order, dim_sf, BasisType::GaussLegendre);
      ParFiniteElementSpace Ws(&vis_mesh, &l2_space_fec);

      L2_FECollection l2_time_fec(order_l2, 1, BasisType::GaussLegendre);
      ParFiniteElementSpace Wt(&time_mesh, &l2_time_fec);

      ParFiniteElementSpace WVs(&vis_mesh, &l2_space_fec, dim_s);

      const int nWs = Ws.GetTrueVSize();
      const int nWt = Wt.GetTrueVSize();

      GridFunction time_nodes(&Wt); // time coordinates
      time_mesh.GetNodes(time_nodes);

      ParGridFunction rho_gf(&Ws);
      ParGridFunction vel_gf(&WVs);
      ParGridFunction mom_gf0(&Ws);

      Vector tdof_vec(nWs);

      string folder = "ParaView/surface/";
      folder += "TestCase"+to_string(testCase);
      folder += "/" + to_string(Vs.GlobalTrueVSize()) + 'x' + to_string(nVt);

      for (int k = 0; k<N;k++){
          ParaViewDataCollection pv(lor ? "SURF_LOR"+to_string(k) : "SURF"+to_string(k), &vis_mesh);
          pv.SetPrefixPath(folder);
          pv.SetHighOrderOutput(!lor);
          pv.SetLevelsOfDetail(vis_order + 1);
          pv.RegisterField("rho", &rho_gf);
          pv.RegisterField("vel", &vel_gf);
	      pv.SetDataFormat(VTKFormat::BINARY32);
          pv.SetCompression(true);

          auto save_time = [&](const double t)
          {
             GetTimeSlice(u_qf, tdof_vec, t, dim_s+k*dim, dim*N, Ws, Wt, lor);
             rho_gf.SetFromTrueDofs(tdof_vec);

             for (int d=0; d<dim_s; d++)
             {
                GetTimeSlice(u_qf, tdof_vec, t, d+k*dim, dim*N, Ws, Wt, lor);
                mom_gf0.SetFromTrueDofs(tdof_vec);
                for (int idx =0; idx< nWs; idx++)
                {
                   double val =  mom_gf0[idx]/rho_gf[idx];
                   vel_gf[idx + d*nWs] = (abs(val) < 10) ? val : 0.0;
                }
             }

             pv.SetTime(t);
             pv.SetCycle(pv.GetCycle() + 1);
             pv.Save();
          };

          for (int k = 0; k < 11; ++k)
          {
              save_time(0.1*k);
          }
      }
   };
   
   // visualization initial data
   auto vis_init = [&](bool lor)
   {
      ParMesh lor_space_mesh;
      if (lor)
      {
         lor_space_mesh = ParMesh::MakeRefined(space_mesh, order_h1, Quadrature1D::GaussLobatto);
      }
      ParMesh &vis_mesh = lor ? lor_space_mesh : space_mesh;
      const int vis_order = lor ? 0 : order_l2;

      L2_FECollection l2_space_fec(vis_order, dim_sf, BasisType::GaussLegendre);
      ParFiniteElementSpace Ws(&vis_mesh, &l2_space_fec);

      L2_FECollection l2_time_fec(order_l2, 1, BasisType::GaussLegendre);
      ParFiniteElementSpace Wt(&time_mesh, &l2_time_fec);

      const int nWs = Ws.GetTrueVSize();
      const int nWt = Wt.GetTrueVSize();

      GridFunction time_nodes(&Wt); // time coordinates
      time_mesh.GetNodes(time_nodes);

      ParGridFunction rho_gf(&Ws);

      Vector tdof_vec(nWs);

      // TODO change this file name..
      string folder = "ParaView/surface/";
    
      folder += "TestCase"+to_string(testCase);
      folder += "/" + to_string(Vs.GlobalTrueVSize()) + 'x' + to_string(nVt);

      for (int k = 0; k<N;k++){
          ParaViewDataCollection pv(lor ? "INIT_SURF_LOR"+to_string(k) : "INIT_SURF"+to_string(k), &vis_mesh);
          pv.SetPrefixPath(folder);
          pv.SetHighOrderOutput(!lor);
          pv.SetLevelsOfDetail(vis_order + 1);
          pv.RegisterField("rho", &rho_gf);
	      pv.SetDataFormat(VTKFormat::BINARY32);
          pv.SetCompression(true);

          auto save_time = [&](const double t)
          {
             GetTimeSlice(u_qf, tdof_vec, t, dim_s+k*dim, dim*N, Ws, Wt, lor);
             rho_gf.SetFromTrueDofs(tdof_vec);
             pv.SetTime(t);
             pv.SetCycle(pv.GetCycle() + 1);
             pv.Save();
          };
          save_time(0.5);
      }
   };
   
   // visualization for terminal density (the barycenter)
   auto vis_bary = [&](bool lor)
   {
      ParMesh lor_space_mesh;
      if (lor)
      {
         lor_space_mesh = ParMesh::MakeRefined(space_mesh, order_h1, Quadrature1D::GaussLobatto);
      }
      ParMesh &vis_mesh = lor ? lor_space_mesh : space_mesh;
      const int vis_order = lor ? 0 : order_l2;

      L2_FECollection l2_space_fec(vis_order, dim_sf, BasisType::GaussLegendre);
      ParFiniteElementSpace Ws(&vis_mesh, &l2_space_fec);

      L2_FECollection l2_time_fec(order_l2, 1, BasisType::GaussLegendre);
      ParFiniteElementSpace Wt(&time_mesh, &l2_time_fec);

      const int nWs = Ws.GetTrueVSize();
      const int nWt = Wt.GetTrueVSize();

      GridFunction time_nodes(&Wt); // time coordinates
      time_mesh.GetNodes(time_nodes);

      ParGridFunction rho_gf(&Ws);

      Vector tdof_vec(nWs);

      string folder = "ParaView/surface/";

      folder += "TestCase"+to_string(testCase);
      folder += "/" + to_string(Vs.GlobalTrueVSize()) + 'x' + to_string(nVt);

      ParaViewDataCollection pv(lor ? "TERM_SURF_LOR": "TERM_SURF", &vis_mesh);
      pv.SetPrefixPath(folder);
      pv.SetHighOrderOutput(!lor);
      pv.SetLevelsOfDetail(vis_order + 1);
      pv.RegisterField("rho", &rho_gf);
	  pv.SetDataFormat(VTKFormat::BINARY32);
      pv.SetCompression(true);

      rho_gf = rhoT_qf; // terminal density
      pv.SetTime(0.5);
      pv.SetCycle(pv.GetCycle() + 1);
      pv.Save();
   };
   
   // store error
   vector<double> errList;
   // TODO 1: SAVE INIT DATA
   vis_init(true);

   for (int it = 0; it < iterALGMax; ++it)
   {
      // Step 1: phi updates
      // initialize delta phi to zero
      sw_a.Start();
      phi_qf = 0.0;
      phiT_qf = 0.0;
      
      for (int k=0; k<N; k++){
        // the (dynamic) linear form for phi
        B.Update(u_qf, s_qf, phi_qf, k, N, sigma_phi);
        // Add the boundary contributions
        {
           Vector B_bdr_slice_0(B, 0, nVs);
           B_bdr.GetColumn(k, B_bdr_0);
           B_bdr_slice_0 += B_bdr_0;

           // Terminal boundary needs to be recomputed
           set_rhs_qf(k);
           b_bdr_T.Assemble();
           b_bdr_T.ParallelAssemble(B_bdr_T);
           Vector B_bdr_slice_T(B, n_time*nVs, nVs);
           B_bdr_slice_T += B_bdr_T;
        }
        
        X0 = 0.0; // initialize to zero
        cg.Mult(B,X0);
        cgIt[k] = cg.GetNumIterations();
        cgErr[k] = cg.GetFinalNorm();
        dX.SetCol(k, X0);
        //dX.Print(cout);

        // Interpolate values at vol & top bdry quadrature points
        phi_qf0.ProjectValue(X0, Vs, Vt);
        face_interp.Project(X0, phiT_qf0, n_time);
        // update phi (incremental)
        for (int i =0; i<nW; i++){
          phi_qf[i*N+k] = phi_qf0[i];
        }

        // update phiT (incremental)
        for (int i=0; i<nWf; i++){
            phiT_qf[i*N+k] = phiT_qf0[i];
        }
      }
      
      // from incremental -> next phi 
      X += dX;
      dX += X;
      ds_qf = 0.0;
      drhoT_qf = 0.0;
      for (int k=0; k<N; k++){
        dX.GetColumn(k, X0);
        phi_qf0.ProjectValue(X0, Vs, Vt);
        dphi_qf0.ProjectGradient(X0, Vs, Vt);
        face_interp.Project(X0, phiT_qf0, n_time);
        // update du_qf: the space-time gradient 
        for (int i =0; i<nW; i++){
          for (int j=0; j<dim; j++){
              du_qf[i*dim*N+dim*k+j] = dphi_qf0[i*dim+j];
          }
        }
        // update ds_qf: the pairwise reaction term
        // FIXME: this is a special example of cylic reaction
        for (int i =0; i<nW; i++){
            ds_qf[i*N+k] += phi_qf0[i];
            ds_qf[i*N+(k+N-1)%N] -= phi_qf0[i];
        }

        // update drhoT_qf : sigma_u * (negative of sum phiT_qf)
        for (int i =0; i<nWf; i++){
            drhoT_qf[i] -= sigma_u * phiT_qf0[i];
        }
      }
      
      sw_a.Stop();
      rt[0] = sw_a.RealTime();
      
      // *qf1 is the incremental
      //u_qf1.Set(-1.0,  u_qf);
      //s_qf1.Set(-1.0,  s_qf);
      rhoT_qf1.Set(-1.0,  rhoT_qf);
      sw_b.Start();
      // STEP 2: nonlinear solver for u_qf
      Vector rho(N), rhobars(N), mbars2(N), sbars(N), sbars2(N), mbars((dim-1)*N); 
      int iter=0;
      for (int j = 0; j < nW; j++)
      {
           // previous data
           mbars2 = 0.0;
           for (int k=0; k<N; k++){
             int idr = j*dim*N+k*dim+dim-1;
             rhobars(k) = sigma_u * du_qf(idr) + u_qf(idr);
             rho(k) = u_qf(idr);// initial guess
             for (int d=0; d<dim-1; d++){
                int idx = j*dim*N + k*dim + d;
                mbars(k*(dim-1)+d) = sigma_u * du_qf(idx) + u_qf(idx);
                mbars2(k) += pow(mbars(k*(dim-1)+d), 2);
             }
           }
         
           // 3 reactions
           for (int k=0; k<N; k++) {
             sbars(k) = sigma_u * ds_qf(j*N+k) + s_qf(j*N+k);
             sbars2(k) = pow(sbars(k),2);
           }
           
           // We use SGS iteration
           for (int k=0; k<N; k++){
             std::uintmax_t it = ptI;
             double mbar2 = mbars2[k];
             double rhobar = rhobars[k];
             auto func = [rho, rhobar, mbar2, sbars2, k](double x) { return 
               F(x, rho, rhobar, mbar2, sbars2, k); };
             auto res = boost::math::tools::brent_find_minima(func, rMin, rMax,
                 double_bits, it);
             rho(k) = res.first;
             if (it > iter) iter = it;
           }

           for (int k=0; k<N; k++){
             for (int d = 0; d < dim; d++)
             {
               if (d<dim-1){
                   auto v1 = V1(rho(k));
                   u_qf(j*dim*N + k*dim + d) = v1/(sigma_u+v1)*mbars(k*(dim-1)+d);
               }else
                   u_qf(j*dim*N + k*dim + d) = rho(k);
             }
           }
           
           // update source
           for (int k=0; k<N;k++){
               auto v2a=V2(rho(k), rho((k+1)%N), alpha);
               s_qf(j*N+k) = v2a/(sigma_u+v2a)*sbars(k);
           }
      }

      // STEP 3: terminal density update (TODO: positivity)
      for (int j=0; j < nWf; j++)
        rhoT_qf[j] = max(rMin, rhoT_qf[j]+drhoT_qf[j]);

      // This is the incremental
      //u_qf1 += u_qf;
      //s_qf1 += s_qf;
      rhoT_qf1 += rhoT_qf;

      sw_b.Stop();
      rt[1] = sw_b.RealTime();

      // compute L1-error for rhoT_qf
      double err = rhoT_qf1.Normlinf();
      MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &iter,  1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      
      // use relative tolerance
      if (it==0)
        tolPDHG *= err;

      if (myid ==0){
          errList.push_back(err);
          if (myid==0){
              cout << std::scientific ;
              cout << "ALG step " << std::setw(3) << it
                 << ",  "<< std::setprecision(3) <<"CG: ";
              for(int k=0;k<N;k++){
                  cout << cgIt[k]<<" "
                       << cgErr[k]<<" ";
              }
              cout << ",  "<< std::setprecision(3) <<"PNT: "<< iter 
                 << ",  "<< std::setprecision(3) <<"err: "<< err
                 << ",  "<<  std::setw(8) << std::setprecision(3) << rt[0] 
                 << ",  "<<  std::setw(8) << std::setprecision(3) << rt[1] 
                 << endl;
          }
      }
 
      // save data
      if ((it+1)%iterPnt==0 && paraview){
          vis(false);
          vis(true);
          vis_bary(true);
          if (myid==0){
             std::ofstream outfile(file);
             if (outfile.is_open()) {
                 for (const auto& err : errList) {
                     outfile << err << std::endl;
                 }
                 outfile.close();
                 std::cout << "Err saved to " << file <<  std::endl;
             } else {
                 std::cerr << "Unable to open file for writing." << std::endl;
                 return 1;
             }
          }
             }
             if (err < tolPDHG)
                break;
   }
   
   return 0;
}


double rhoInit(const Vector &x, double t, int k)
{
  if (testCase==1){// gauss translations in 3D
    double xc = 1.0;
    double yc = 0.0;
    double zc = k-0.5;
    double r2 = pow(x(0)-xc,2)+pow(x(1)-yc,2) + pow(x(2)-zc, 2);
    double val = exp(-20*r2);
    return val;
  }else if (testCase==2){// 3 gaussian translation
    double theta = 1.5*M_PI*(k-1)/N;
    // Gaussian centers
    double xc = x(0) - cos(theta);
    double yc = x(1) - sin(theta);
    double zc = x(2) - 0.0;
    double r2;
    if (k==0)
       r2 = xc*xc + yc* yc + 2.0*zc*zc + 2 * xc * zc;
    if (k==1)
       r2 = xc*xc + 2.0*yc* yc  + zc*zc;
    if (k==2)
       r2 = 2.0*xc*xc + yc* yc + zc*zc - 2 * xc * zc;

    double val = exp(-20*r2);
    return val;
  }else if (testCase==3){// armadilo
    double val = 1e-6; // background val.
    if (k==0)// left foot
      if (x(1) < -0.12 && x(0) >0) val = 1.0;
    if (k==1)// right foot
      if (x(2) > 0.64 ) val = 1.0;
    if (k==2)// left hand
      if (x(1) > 0.80 && x(0) > 0.30) val = 1.0;
    if (k==3)// right hand
      if (x(1) > 0.78 && x(0) < -0.40) val = 1.0;
    if (k==4)// ears
      if (x(1) > 0.92 && x(0) < 0.30 && x(0) > -0.30 ) val = 1.0;

    return val;
  }else if (testCase==4){// dinosaur
    double val = 1e-6; // background val.
    if (k==0){// head
        if (x(1)>0.55) val = 1.0; 
    }
    if (k==1){// tail
        if (x(1)< - 0.45) val = 1.0; 
    }
    if (k==2){// back left toe
        if (x(1) < -0.40 && x(2) < 0.01) val = 1.0; 
    }
    if (k==3){// front right toe
        if (x(0) > 0.40) val = 1.0; 
    }

    return val;
  }else
    return 0.25; // default choice
}

double V1(double rho){ return rho; }

double V2(double rho1, double rho2, double c){
   double bot = log(rho1)-log(rho2);
   double top = (rho1-rho2);
   double val;
   if (abs(bot)<1e-8) val = rho2;
   else
    val = top/bot;
   return c*val;
}

// interaction potential (serves as a regularizer)
double E(double rho){ return rho*log(rho); }

double F(double rho, Vector rhos,  double rhobar, double mbar2, 
    Vector sbars2, int k){
   int idX0=(k+1)%N, idX1 = (k+N-1)%N, idY0 = k, idY1 = (k+N-1)%N;
   auto e =   E(rho);
   auto v1 =   V1(rho);
   auto v2A =    V2(rho, rhos(idX0), alpha);
   auto v2B =    V2(rhos(idX1), rho, alpha);

   return 0.5*pow(rho - rhobar,2)/sigma_u 
       + 0.5*mbar2/(sigma_u + v1)
       + 0.5*sbars2(idY0)/(sigma_u + v2A)
       + 0.5*sbars2(idY1)/(sigma_u + v2B)
       + beta*e;
}

// NOT WORKING on surface mesh yet
double L1Norm(QuadratureFunction qf)
{
   const auto Qs = qf.GetSpace();
   const Vector &ws = Qs->GetWeights();

   const int nQs = Qs->GetSize();

   double integ = 0.0;
   for (int is = 0; is < nQs; ++is)
   {
      const double val = qf[is];
      integ += ws[is]*std::abs(val);
   }

   if (auto *pmesh = dynamic_cast<const ParMesh*>(Qs->GetMesh()))
   {
      MPI_Comm comm = pmesh->GetComm();
      MPI_Allreduce(MPI_IN_PLACE, &integ, 1, MPI_DOUBLE, MPI_SUM, comm);
   }
   return integ;
}
