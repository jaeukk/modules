/** \file experiment.cpp
 * \brief Contact-energy reproducibility controls using the paired experiment harness.
 * Compile CC_TEST_N/SPECIES/THREADS/MODES to select the matched benchmark.
 */
#define main paired_experiment_main
#include "../five_component/experiment.cpp"
#undef main

/** \brief Check independent forces, finite differences, and repeated fixed-point values. */
void probe() {
    auto c=config(20260928);auto m=modes(Species==1?0.3:0.1);
    auto ref=reference(c,m);
    std::cout<<std::setprecision(17)<<"{\"N\":"<<N<<",\"species\":"<<Species<<",\"threads\":"<<Threads<<",\"backends\":[";
    bool first=true;
    for(auto backend:{Backend::Direct,Backend::FinufftType1}){
        std::unique_ptr<MultiStealthyPotential> p;
        {Quiet quiet;p=std::make_unique<MultiStealthyPotential>(2,1);setup(*p,c,m,backend);}
        std::vector<GeometryVector> initial,current;p->AllForce(initial);
        const double referenceError=forceDifference(initial,ref.force);
        check(referenceError<1e-8,"independent force check failed");
        double lo=1e100,hi=-1e100,maxForceChange=0,maxFD=0;
        for(int j=0;j<100;++j){
            p->SetConfiguration(c);double e=p->Energy();lo=std::min(lo,e);hi=std::max(hi,e);
            check(std::abs(e-ref.cc-ref.contact)<1e-8,"independent energy check failed");
            p->AllForce(current);maxForceChange=std::max(maxForceChange,forceDifference(initial,current));
            if(backend==Backend::FinufftType1)for(auto& s:p->species_constraints)
                check(s.CCPotential->GetLastRhoBackendUsed()==Backend::FinufftType1,"probe fallback");
        }
        for(size_t i:{size_t(0),N/3,N-1})for(int d=0;d<2;++d){
            double e[2],h=1e-5;
            for(int sign=0;sign<2;++sign){auto moved=c;auto s=c.GetRelativeCoordinates(i);
                s.x[d]+=(sign?-h:h)/L;moved.MoveParticle(i,s);p->SetConfiguration(moved);e[sign]=p->Energy();}
            maxFD=std::max(maxFD,std::abs(-(e[0]-e[1])/(2*h)-initial[i].x[d]));
        }
        check(maxFD<3e-6,"finite difference check failed");
        if(!first)std::cout<<',';first=false;
        std::cout<<"{\"backend\":\""<<(backend==Backend::Direct?"direct":"finufft")<<"\",\"energy_range\":"<<hi-lo
                 <<",\"force_range\":"<<maxForceChange<<",\"reference_force_error\":"<<referenceError<<",\"finite_difference_error\":"<<maxFD<<'}';
    }
    std::cout<<"],\"passed\":true}\n";
}

/** \brief Run either the nonoptimizing probe or an ordinary saved optimization. */
int main(int argc,char** argv)try{
    if(argc==2&&std::string(argv[1])=="probe"){probe();return 0;}
    return paired_experiment_main(argc,argv);
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
