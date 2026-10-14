// Five-component gate and finite random-start experiment; outputs are reproducible artifacts.
#include "MultiStealthyPotential.h"
#include "StructureOptimization.h"
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <random>
#include <sstream>

using Backend=ShiftedCCPotential::RhoBackend;
using Clock=std::chrono::steady_clock;
#ifndef CC_TEST_N
#define CC_TEST_N 2000
#endif
#ifndef CC_TEST_SPECIES
#define CC_TEST_SPECIES 5
#endif
#ifndef CC_TEST_THREADS
#define CC_TEST_THREADS 4
#endif
#ifndef CC_TEST_MODES
#define CC_TEST_MODES 0
#endif
constexpr size_t N=CC_TEST_N, Species=CC_TEST_SPECIES, Ns=N/Species;
constexpr int Threads=CC_TEST_THREADS;
static_assert(N % Species == 0, "equal component sizes required");
const double L=std::sqrt(double(N)), Radius=std::sqrt(0.2/::pi), Diameter=2*Radius;
size_t Verbosity=0;
double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
struct Pass{int status,count;double energy;};
std::vector<Pass> passes;
Clock::time_point optimizationStart;
extern "C" nlopt_result __real_nlopt_optimize(nlopt_opt,double*,double*);
extern "C" nlopt_result __wrap_nlopt_optimize(nlopt_opt opt,double* x,double* f){
    double remaining=600-elapsed(optimizationStart);
    if(remaining<=0){*f=passes.back().energy;passes.push_back({6,0,*f});return NLOPT_MAXTIME_REACHED;}
    nlopt_set_maxtime(opt,remaining);
    auto status=__real_nlopt_optimize(opt,x,f);
    passes.push_back({int(status),nlopt_get_numevals(opt),*f});return status;
}
struct Quiet{std::ostringstream sink;std::streambuf* old;
    Quiet():old(std::cout.rdbuf(sink.rdbuf())){} ~Quiet(){std::cout.rdbuf(old);}};
struct Modes{std::vector<std::array<int,2>> indices;std::vector<GeometryVector> k;int shell;};
Modes modes(double target){
    std::map<int,std::vector<std::array<int,2>>> shells;
    for(int y=-25;y<=25;++y)for(int x=-25;x<=25;++x)
        if(y>0||(y==0&&x>0))shells[x*x+y*y].push_back({x,y});
    Modes best;double error=1e30;std::vector<std::array<int,2>> accumulated;
    for(const auto& shell:shells){
        accumulated.insert(accumulated.end(),shell.second.begin(),shell.second.end());
        double delta=std::abs(double(accumulated.size())/(2*(Ns-1))-target);
        if(delta<error){error=delta;best.indices=accumulated;best.shell=shell.first;}
    }
    if(CC_TEST_MODES>0){
        best.indices.clear();
        for(const auto& shell:shells)for(auto index:shell.second)
            if(best.indices.size()<CC_TEST_MODES){best.indices.push_back(index);best.shell=shell.first;}
        check(best.indices.size()==CC_TEST_MODES,"insufficient test modes");
    }
    for(auto index:best.indices){GeometryVector k(2);k.x[0]=2*::pi*index[0]/L;k.x[1]=2*::pi*index[1]/L;best.k.push_back(k);}
    return best;
}
Configuration config(unsigned long long seed,bool lattice=false){
    GeometryVector basis[2]={GeometryVector(2),GeometryVector(2)};basis[0].x[0]=L;basis[1].x[1]=L;
    Configuration c(2,basis,1.0);std::mt19937_64 rng(seed);std::uniform_real_distribution<double> u(0,1);
    for(size_t i=0;i<N;++i){GeometryVector s(2);
        if(lattice){size_t cell=i/Species,sp=i%Species;s.x[0]=(cell%20+double(sp)/5)/20;s.x[1]=(cell/20+double(sp)/5)/20;}
        else{s.x[0]=u(rng);s.x[1]=u(rng);}c.Insert("a",s);}
    return c;
}
void save(const Configuration& c,const std::filesystem::path& path){
    std::ofstream out(path);check(bool(out),"cannot save coordinates");out<<std::setprecision(17)<<"id\tspecies\tsx\tsy\tradius\n";
    for(size_t i=0;i<N;++i){auto s=c.GetRelativeCoordinates(i);out<<i<<'\t'<<i%Species<<'\t'<<s.x[0]<<'\t'<<s.x[1]<<'\t'<<Radius<<'\n';}
}
void setup(MultiStealthyPotential& p,const Configuration& c,const Modes& m,Backend backend,bool reverse=false){
    p.ParallelNumThread=Threads;p.SetConfiguration(c);p.SetNumSpecies(Species);p.SetSecureRadii(std::vector<double>(N,Radius));
    for(size_t j=0;j<N;++j){size_t i=reverse?N-1-j:j;p.SetSpeciesOfParticle(i,i%Species);}
    for(size_t s=0;s<Species;++s){for(const auto& k:m.k)p.AddConstraint(s,k,{1.0});p.species_constraints[s].CCPotential->SetRhoBackend(backend);}
    p.SetConfiguration(c);
}
struct Counted:ShiftedCCPotential{
    size_t count=0,fast=0;double seconds=0;
    explicit Counted(const ShiftedCCPotential& p):ShiftedCCPotential(p){}
    // Reads particles only through the base class, so MultiStealthyPotential may bind its indexed view.
    bool SupportsIndexedView()const override{return true;}
    void GetRho()override{if(RhoValid)return;auto t=Clock::now();ShiftedCCPotential::GetRho();seconds+=elapsed(t);++count;fast+=GetLastRhoBackendUsed()==Backend::FinufftType1;}
};
struct Reference{
    double cc=0,contact=0,maxOverlap=0,maxForce=0,maxS=0;
    std::array<double,Species> speciesS{};
    std::vector<GeometryVector> force;
    bool ground()const{return cc+contact<=1e-11&&maxS<=1e-10&&maxForce<=1e-7&&maxOverlap<=1e-6;}
};
Reference reference(const Configuration& c,const Modes& m,const std::filesystem::path& file={}){
    Reference r;r.force.assign(N,GeometryVector(2));std::ofstream out;
    if(!file.empty()){out.open(file);check(bool(out),"cannot save structure factors");out<<std::setprecision(17)<<"species\tmx\tmy\tS\n";}
    for(size_t s=0;s<Species;++s)for(size_t index=0;index<m.k.size();++index){
        const auto& k=m.k[index];std::array<double,Ns> co,si;double re=0,im=0;
        for(size_t j=0;j<Ns;++j){double phase=k.Dot(c.GetCartesianCoordinates(Species*j+s));co[j]=std::cos(phase);si[j]=std::sin(phase);re+=co[j];im+=si[j];}
        double ss=(re*re+im*im)/Ns;r.speciesS[s]=std::max(r.speciesS[s],ss);r.cc+=(re*re+im*im)/(Species*L*L);
        if(out.is_open())out<<s<<'\t'<<m.indices[index][0]<<'\t'<<m.indices[index][1]<<'\t'<<ss<<'\n';
        for(size_t j=0;j<Ns;++j)r.force[Species*j+s].AddFrom((2.0/(Species*L*L)*(re*si[j]-im*co[j]))*k);
    }
    for(size_t i=0;i<N;++i)for(size_t j=i+1;j<N;++j){
        auto delta=c.GetRelativeCoordinates(j)-c.GetRelativeCoordinates(i);
        for(int d=0;d<2;++d)delta.x[d]=L*(delta.x[d]-std::round(delta.x[d]));
        double rr=delta.Modulus2();if(rr>=Diameter*Diameter)continue;
        double distance=std::sqrt(rr);check(distance>1e-14,"coincident particles in reference");
        double overlap=1-distance/Diameter;r.contact+=overlap*overlap;r.maxOverlap=std::max(r.maxOverlap,overlap);
        auto f=(-2.0/Diameter*(1/distance-1/Diameter))*delta;r.force[i].AddFrom(f);r.force[j].MinusFrom(f);
    }
    for(auto f:r.force)for(int d=0;d<2;++d){check(std::isfinite(f.x[d]),"nonfinite reference force");r.maxForce=std::max(r.maxForce,std::abs(f.x[d]));}
    r.maxS=*std::max_element(r.speciesS.begin(),r.speciesS.end());check(std::isfinite(r.cc+r.contact),"nonfinite energy");return r;
}
void print(const Reference& r){
    std::cout<<"{\"cc_energy\":"<<r.cc<<",\"contact_energy\":"<<r.contact<<",\"total_energy\":"<<r.cc+r.contact
      <<",\"max_relative_overlap\":"<<r.maxOverlap<<",\"max_force_component\":"<<r.maxForce<<",\"max_S\":"<<r.maxS<<",\"species_max_S\":[";
    for(size_t s=0;s<Species;++s){if(s)std::cout<<',';std::cout<<r.speciesS[s];}
    std::cout<<"],\"ground_state\":"<<(r.ground()?"true":"false")<<'}';
}
double forceDifference(const std::vector<GeometryVector>& a,const std::vector<GeometryVector>& b){
    double result=0;for(size_t i=0;i<N;++i)for(int d=0;d<2;++d)result=std::max(result,std::abs(a[i].x[d]-b[i].x[d]));return result;
}
void validate(const Modes& m){
    Configuration c=config(20260928);double maxFD=0,backendForce=0,backendRho=0,refForce=0,energyError=0;
    {Quiet quiet;MultiStealthyPotential direct(2,1),fast(2,1),reverse(2,1);
    setup(direct,c,m,Backend::Direct);setup(fast,c,m,Backend::FinufftType1);setup(reverse,c,m,Backend::Direct,true);
    auto ref=reference(c,m);std::vector<GeometryVector> df,ff,rf;
    // Cold force evaluation must fill shared caches before starting its parallel loop.
    direct.AllForce(df);fast.AllForce(ff);reverse.AllForce(rf);
    check(forceDifference(df,rf)<1e-10,"species assignment order changes force");
    backendForce=forceDifference(df,ff);refForce=forceDifference(df,ref.force);
    energyError=std::abs(fast.Energy()-direct.Energy());
    check(std::abs(direct.Energy()-ref.cc-ref.contact)<1e-8,"energy reference mismatch");
    check(backendForce<1e-8&&refForce<1e-8&&energyError<1e-8,"backend/reference comparison failed");
    for(size_t s=0;s<Species;++s){auto& d=*direct.species_constraints[s].CCPotential;auto& f=*fast.species_constraints[s].CCPotential;
        check(f.GetLastRhoBackendUsed()==Backend::FinufftType1,"validation fell back");
        for(size_t k=0;k<m.k.size();++k)backendRho=std::max({backendRho,std::abs(d.RhoReal[k]-f.RhoReal[k]),std::abs(d.RhoImag[k]-f.RhoImag[k])});}
    check(backendRho<1e-8,"density backend mismatch");
    // Reassign out and back, changing insertion order without changing final labels.
    reverse.SetSpeciesOfParticle(0,1);reverse.SetSpeciesOfParticle(0,0);reverse.AllForce(rf);
    check(forceDifference(df,rf)<1e-10,"species reassignment changes force");
    for(auto* p:{&direct,&fast})for(size_t i:{size_t(0),size_t(17),size_t(333),size_t(999),size_t(1500),size_t(1999)})for(int d=0;d<2;++d){
        double e[2],h=1e-5;for(int sign=0;sign<2;++sign){Configuration moved(c);auto s=c.GetRelativeCoordinates(i);s.x[d]+=(sign?-h:h)/L;moved.MoveParticle(i,s);p->SetConfiguration(moved);e[sign]=p->Energy();}
        double err=std::abs(-(e[0]-e[1])/(2*h)-df[i].x[d]);maxFD=std::max(maxFD,err);
        check(err<3e-6*(1+std::abs(df[i].x[d])),"combined force finite difference failed");}
    // Known exact zero-energy fixture, five translated 20x20 lattices.
    Configuration lattice=config(0,true);
    auto lr=reference(lattice,m);check(lr.ground(),"constructed reference is not ground state");
    for(auto* p:{&direct,&fast}){p->SetConfiguration(lattice);check(p->Energy()<1e-18,"constructed backend energy failed");p->AllForce(ff);check(forceDifference(ff,lr.force)<1e-8,"constructed force failed");}
    // Radius setter must refresh a previously too-small neighbor cutoff.
    GeometryVector basis[2]={GeometryVector(2),GeometryVector(2)};basis[0].x[0]=4;basis[1].x[1]=4;
    Configuration two(2,basis,0.5);GeometryVector a(2),b(2);b.x[0]=0.3;two.Insert("a",a);two.Insert("a",b);
    MultiStealthyPotential contact(2,1);contact.SetConfiguration(two);contact.SetNumSpecies(1);contact.SetSpeciesOfParticle(0,0);contact.SetSpeciesOfParticle(1,0);contact.SetSecureRadii({0.7,0.7});
    check(std::abs(contact.Energy()-std::pow(1-1.2/1.4,2))<1e-12,"radius cutoff stale");}
    std::cout<<"{\"passed\":true,\"N\":"<<N<<",\"threads\":4,\"mode_pairs\":"<<m.k.size()<<",\"actual_chi\":"<<double(m.k.size())/(2*(Ns-1))
      <<",\"max_density_difference\":"<<backendRho<<",\"energy_difference\":"<<energyError<<",\"force_backend_difference\":"<<backendForce
      <<",\"force_reference_difference\":"<<refForce<<",\"finite_difference_max_error\":"<<maxFD<<",\"constructed_checks\":[";
    bool first=true;Configuration lattice=config(0,true);
    for(double target:{0.1,0.2,0.3,0.4,0.45,0.5,0.55}){auto mm=modes(target);auto r=reference(lattice,mm);check(r.ground(),"constructed sweep reference failed");
        {Quiet quiet;MultiStealthyPotential fast(2,1);setup(fast,lattice,mm,Backend::FinufftType1);check(fast.Energy()<1e-18,"constructed FINUFFT sweep check failed");}
        if(!first)std::cout<<',';first=false;std::cout<<"{\"target_chi\":"<<target<<",\"mode_pairs\":"<<mm.k.size()<<",\"actual_chi\":"<<double(mm.k.size())/(2*(Ns-1))<<",\"reference\":";print(r);std::cout<<'}';}
    std::cout<<"]}\n";
}
int main(int argc,char** argv)try{
    std::cout<<std::setprecision(17);TimeLimit=std::numeric_limits<time_t>::max();StructureOptimization_SaveConfigurationInterval=1000000000;
    if(argc==2&&std::string(argv[1])=="validate"){validate(modes(0.1));return 0;}
    check(argc==5,"usage: experiment OUTDIR SEED CHI direct|finufft");std::filesystem::path out=argv[1];
    check(!std::filesystem::exists(out),"refusing to overwrite an experiment");std::filesystem::create_directories(out);
    auto seed=std::stoull(argv[2]);double target=std::stod(argv[3]);std::string backend=argv[4];check(backend=="direct"||backend=="finufft","invalid backend");
    auto m=modes(target);Configuration c=config(seed);save(c,out/"initial.tsv");
    std::ofstream mf(out/"modes.tsv");for(auto index:m.indices)mf<<index[0]<<'\t'<<index[1]<<'\n';mf.close();
    std::unique_ptr<MultiStealthyPotential> p;
    {Quiet quiet;p=std::make_unique<MultiStealthyPotential>(2,1);setup(*p,c,m,backend=="direct"?Backend::Direct:Backend::FinufftType1);
        for(auto& s:p->species_constraints)s.CCPotential=std::make_unique<Counted>(*s.CCPotential);
        p->SetConfiguration(c);p->Energy();}
    for(auto& s:p->species_constraints){auto& q=dynamic_cast<Counted&>(*s.CCPotential);q.count=q.fast=0;q.seconds=0;}
    optimizationStart=Clock::now();RelaxStructure_NLOPT(c,*p,0,0,0,2000);double wall=elapsed(optimizationStart);
    save(c,out/"final.tsv");auto ref=reference(c,m,out/"direct_structure.tsv");
    size_t rhoCalls=0,fastCalls=0;double rhoTime=0;for(auto& s:p->species_constraints){auto& q=dynamic_cast<Counted&>(*s.CCPotential);rhoCalls+=q.count;fastCalls+=q.fast;rhoTime+=q.seconds;}
    check(backend=="direct"||rhoCalls==fastCalls,"FINUFFT fallback during optimization");
    {Quiet quiet;MultiStealthyPotential vanilla(2,1);setup(vanilla,c,m,Backend::Direct);std::vector<GeometryVector> f;double e=vanilla.Energy();vanilla.AllForce(f);
        check(std::abs(e-ref.cc-ref.contact)<1e-8*(1+std::abs(e)),"final energy independent check failed");check(forceDifference(f,ref.force)<1e-8,"final force independent check failed");}
    std::cout<<"{\"seed\":"<<seed<<",\"backend\":\""<<backend<<"\",\"N\":"<<N<<",\"threads\":"<<Threads<<",\"packing_fraction\":0.2,\"target_chi\":"<<target
      <<",\"actual_chi\":"<<double(m.k.size())/(2*(Ns-1))<<",\"mode_pairs\":"<<m.k.size()<<",\"shell_squared\":"<<m.shell<<",\"optimization_seconds\":"<<wall
      <<",\"rho_seconds\":"<<rhoTime<<",\"rho_calls\":"<<rhoCalls<<",\"finufft_calls\":"<<fastCalls<<",\"passes\":[";
    for(size_t i=0;i<passes.size();++i){if(i)std::cout<<',';auto a=passes[i];std::cout<<"{\"status\":"<<a.status<<",\"evaluations\":"<<a.count<<",\"energy\":"<<a.energy<<'}';}
    std::cout<<"],\"direct_reference\":";print(ref);std::cout<<"}\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
