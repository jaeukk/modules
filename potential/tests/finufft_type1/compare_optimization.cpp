// Paired backend experiment using the repository's unmodified NLopt optimizer.
#include "CollectiveCoordinatePotential.h"
#include "StructureOptimization.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>

using Clock = std::chrono::steady_clock;
using Backend = ShiftedCCPotential::RhoBackend;
size_t Verbosity = 0; // Application-owned global required by StructureOptimization.cpp.
constexpr size_t N = 500, M = 299, MaxEvaluationsPerPass = 2000;
constexpr double ELimit = 1e-11, SLimit = 1e-11, FLimit = 1e-7;
double seconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }
void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }

struct Pass { int status, evaluations; double energy; };
std::vector<Pass> passes;
// Link-time observation only: no optimizer settings or return values are changed.
extern "C" nlopt_result __real_nlopt_optimize(nlopt_opt, double*, double*);
extern "C" nlopt_result __wrap_nlopt_optimize(nlopt_opt opt, double* x, double* f) {
    nlopt_result status = __real_nlopt_optimize(opt, x, f);
    passes.push_back({int(status), nlopt_get_numevals(opt), *f});
    return status;
}

struct CountedPotential : ShiftedCCPotential {
    size_t rhoCalls=0, fastCalls=0, forceCalls=0, energyCalls=0;
    double rhoSeconds=0, forceSeconds=0;
    explicit CountedPotential(const ShiftedCCPotential& p) : ShiftedCCPotential(p) {}
    void GetRho() override {
        if (RhoValid) return;
        auto t=Clock::now();
        ShiftedCCPotential::GetRho();
        rhoSeconds+=seconds(t); ++rhoCalls;
        fastCalls += GetLastRhoBackendUsed()==Backend::FinufftType1;
    }
    double Energy() override { ++energyCalls; return ShiftedCCPotential::Energy(); }
    void AllForce(std::vector<GeometryVector>& f) override {
        // The objective calls Energy first; require cached rho to keep timing disjoint.
        require(RhoValid, "optimizer called forces without cached density");
        auto t=Clock::now(); ShiftedCCPotential::AllForce(f);
        forceSeconds+=seconds(t); ++forceCalls;
    }
};

Configuration emptyConfig() {
    double L=std::sqrt(double(N));
    GeometryVector basis[2]={GeometryVector(2),GeometryVector(2)};
    basis[0].x[0]=L; basis[1].x[1]=L;
    return Configuration(2,basis,10*L);
}
void save(const Configuration& c, const std::filesystem::path& file) {
    std::ofstream out(file); require(bool(out), "cannot save positions");
    out<<std::setprecision(17);
    for(size_t j=0;j<N;++j) {
        const auto& s=c.GetRelativeCoordinates(j);
        out<<s.x[0]<<'\t'<<s.x[1]<<'\n';
    }
    require(bool(out), "position write failed");
}
Configuration load(const std::filesystem::path& file) {
    Configuration c=emptyConfig(); std::ifstream in(file);
    for(size_t j=0;j<N;++j) {
        GeometryVector s(2); require(bool(in>>s.x[0]>>s.x[1]), "invalid initial positions");
        c.Insert("A",s);
    }
    return c;
}
ShiftedCCPotential potential(const Configuration& c, const std::filesystem::path& file) {
    std::vector<std::array<int,2>> modes;
    for(int y=-16;y<=16;++y) for(int x=-16;x<=16;++x)
        if(x*x+y*y<=256 && (y>0 || (y==0 && x>0))) modes.push_back({x,y});
    std::sort(modes.begin(),modes.end(),[](const auto& a,const auto& b){
        int aa=a[0]*a[0]+a[1]*a[1],bb=b[0]*b[0]+b[1]*b[1];
        if(aa!=bb) return aa<bb;
        return a[1]!=b[1]?a[1]<b[1]:a[0]<b[0];
    });
    modes.resize(M);
    ShiftedCCPotential p(2); p.ParallelNumThread=2;
    std::ofstream out(file); require(bool(out),"cannot save mode list");
    for(auto m:modes) {
        GeometryVector k=m[0]*c.GetReciprocalBasisVector(0)+m[1]*c.GetReciprocalBasisVector(1);
        p.AddConstraint(k,{1.0}); out<<m[0]<<'\t'<<m[1]<<'\n';
    }
    return p;
}

struct Validation { double energy=0,maxS=0,meanS=0,maxF=0,rmsF=0,rhoError=0,energyError=0,forceError=0; bool ground=false; };
Validation validate(const Configuration& c,const ShiftedCCPotential& proto,const std::filesystem::path& file) {
    ShiftedCCPotential d(proto),f(proto);
    d.SetRhoBackend(Backend::Direct); f.SetRhoBackend(Backend::FinufftType1);
    d.SetConfiguration(c); f.SetConfiguration(c);
    Validation v; v.energy=d.Energy(); v.energyError=std::abs(v.energy-f.Energy());
    require(f.GetLastRhoBackendUsed()==Backend::FinufftType1,"unexpected FINUFFT fallback");
    std::vector<GeometryVector> df,ff; d.AllForce(df); f.AllForce(ff);
    std::ofstream out(file); require(bool(out),"cannot save direct validation");
    out<<std::setprecision(17)<<"mode\trho_real\trho_imag\tS\n";
    for(size_t k=0;k<M;++k) {
        double sr=(d.RhoReal[k]*d.RhoReal[k]+d.RhoImag[k]*d.RhoImag[k])/N;
        require(std::isfinite(sr),"nonfinite structure factor");
        v.maxS=std::max(v.maxS,sr); v.meanS+=sr/M;
        v.rhoError=std::max({v.rhoError,std::abs(d.RhoReal[k]-f.RhoReal[k]),std::abs(d.RhoImag[k]-f.RhoImag[k])});
        out<<k<<'\t'<<d.RhoReal[k]<<'\t'<<d.RhoImag[k]<<'\t'<<sr<<'\n';
    }
    for(size_t j=0;j<N;++j) for(int k=0;k<2;++k) {
        require(std::isfinite(df[j].x[k]) && std::isfinite(ff[j].x[k]),"nonfinite force");
        v.maxF=std::max(v.maxF,std::abs(df[j].x[k])); v.rmsF+=df[j].x[k]*df[j].x[k]/(2*N);
        v.forceError=std::max(v.forceError,std::abs(df[j].x[k]-ff[j].x[k]));
    }
    v.rmsF=std::sqrt(v.rmsF);
    require(std::isfinite(v.energy) && std::isfinite(v.energyError),"nonfinite energy");
    // Density one and unit mode weights imply E = sum S for this experiment.
    require(std::abs(v.energy-v.meanS*M)<=1e-10*(1+std::abs(v.energy)),"energy/structure-factor inconsistency");
    v.ground=v.energy<=ELimit && v.maxS<=SLimit && v.maxF<=FLimit;
    return v;
}
void print(const Validation& v) {
    std::cout<<"{\"energy\":"<<v.energy<<",\"max_S\":"<<v.maxS<<",\"mean_S\":"<<v.meanS
      <<",\"max_force_component\":"<<v.maxF<<",\"rms_force_component\":"<<v.rmsF
      <<",\"rho_backend_error\":"<<v.rhoError<<",\"energy_backend_error\":"<<v.energyError
      <<",\"force_backend_error\":"<<v.forceError<<",\"ground_state\":"<<(v.ground?"true":"false")<<"}";
}

struct Run { Configuration config; Validation validation; double wall; CountedPotential measured; std::vector<Pass> steps;
    Run(const Configuration& c,const ShiftedCCPotential& p):config(c),wall(0),measured(p){} };
Run optimize(const std::filesystem::path& dir,const ShiftedCCPotential& p,Backend b,int repeat) {
    // Both backends reread the same saved double-precision initial configuration.
    Run r(load(dir/"initial.tsv"),p); r.measured.SetRhoBackend(b);
    r.measured.SetConfiguration(r.config); r.measured.GetRho(); // untimed warm-up
    r.measured.rhoCalls=r.measured.fastCalls=0; r.measured.rhoSeconds=0;
    passes.clear();
    auto start=Clock::now();
    RelaxStructure_NLOPT(r.config,r.measured,0.0,0,0.0,MaxEvaluationsPerPass);
    r.wall=seconds(start); r.steps=passes;
    std::string name=(b==Backend::Direct?"direct_":"finufft_")+std::to_string(repeat);
    save(r.config,dir/(name+".tsv"));
    r.validation=validate(r.config,p,dir/(name+"_structure.tsv"));
    if(b==Backend::FinufftType1) require(r.measured.rhoCalls==r.measured.fastCalls,"optimization silently fell back");
    return r;
}
void print(const Run& r) {
    std::cout<<"{\"wall_seconds\":"<<r.wall<<",\"rho_seconds\":"<<r.measured.rhoSeconds
      <<",\"force_seconds\":"<<r.measured.forceSeconds<<",\"rho_evaluations\":"<<r.measured.rhoCalls
      <<",\"finufft_evaluations\":"<<r.measured.fastCalls<<",\"energy_calls\":"<<r.measured.energyCalls
      <<",\"gradient_calls\":"<<r.measured.forceCalls<<",\"passes\":[";
    for(size_t i=0;i<r.steps.size();++i) { if(i)std::cout<<','; auto p=r.steps[i];
        std::cout<<"{\"status\":"<<p.status<<",\"evaluations\":"<<p.evaluations<<",\"energy\":"<<p.energy<<'}'; }
    std::cout<<"],\"direct_validation\":"; print(r.validation); std::cout<<'}';
}

double wrap(double x) {return x-std::floor(x+0.5);}
void displacement(const Configuration& a,const Configuration& b) {
    // Exact labeled least-squares translation on each circle: enumerate all unwrap cuts.
    std::array<double,2> shift{}; double raw2=0,aligned2=0,maxDistance=0;
    for(int d=0;d<2;++d) {
        std::vector<double> delta; double sum=0,best=std::numeric_limits<double>::infinity();
        for(size_t i=0;i<N;++i) {double v=wrap(b.GetRelativeCoordinates(i).x[d]-a.GetRelativeCoordinates(i).x[d]);
            delta.push_back(v); sum+=v; raw2+=v*v;}
        std::sort(delta.begin(),delta.end());
        for(size_t cut=0;cut<N;++cut) {
            double candidate=(sum+cut)/N,score=0;
            for(double v:delta) {double dv=wrap(v-candidate);score+=dv*dv;}
            if(score<best) {best=score;shift[d]=wrap(candidate);}
        }
    }
    for(size_t i=0;i<N;++i) {double norm=0;for(int d=0;d<2;++d) {
        double v=wrap(b.GetRelativeCoordinates(i).x[d]-a.GetRelativeCoordinates(i).x[d]-shift[d]);norm+=v*v;}
        aligned2+=norm;maxDistance=std::max(maxDistance,std::sqrt(norm*N));}
    // L=sqrt(N); convert fractional displacements to Cartesian distances.
    std::cout<<"{\"raw_periodic_rms\":"<<std::sqrt(raw2)<<",\"aligned_rms\":"<<std::sqrt(aligned2)
      <<",\"aligned_max\":"<<maxDistance<<",\"translation_cartesian\":["<<shift[0]*std::sqrt(N)<<','<<shift[1]*std::sqrt(N)<<"]}";
}

void microbench(const Configuration& c,const ShiftedCCPotential& proto) {
    ShiftedCCPotential p[2]={proto,proto}; p[1].SetRhoBackend(Backend::FinufftType1);
    for(auto& q:p) q.SetConfiguration(c);
    std::vector<GeometryVector> forces;
    std::vector<double> samples[2][3];
    for(int iteration=-5;iteration<30;++iteration) for(int j=0;j<2;++j) {
        int b=(iteration+6+j)%2; auto& q=p[b];
        q.RhoValid=false; auto t=Clock::now();q.GetRho();double rho=seconds(t);
        t=Clock::now();q.AllForce(forces);double force=seconds(t);
        t=Clock::now();q.SetConfiguration(c);double energy=q.Energy();q.AllForce(forces);double full=seconds(t);
        require(std::isfinite(energy),"nonfinite benchmark energy");
        require(q.GetLastRhoBackendUsed()==(b?Backend::FinufftType1:Backend::Direct),"benchmark backend mismatch");
        if(iteration>=0) {samples[b][0].push_back(rho);samples[b][1].push_back(force);samples[b][2].push_back(full);}
    }
    std::cout<<'{';
    for(int b=0;b<2;++b) {if(b)std::cout<<',';std::cout<<'"'<<(b?"finufft":"direct")<<"\":{";
        for(int kind=0;kind<3;++kind) {if(kind)std::cout<<',';std::cout<<'"'<<(kind==0?"rho_seconds":kind==1?"cached_force_seconds":"energy_force_seconds")<<"\":[";
            for(size_t i=0;i<samples[b][kind].size();++i){if(i)std::cout<<',';std::cout<<samples[b][kind][i];}std::cout<<']';}
        std::cout<<'}';}
    std::cout<<'}';
}

int main(int argc,char** argv) try {
    require(argc==4,"usage: compare_optimization OUTDIR SEED REPETITIONS");
    const std::filesystem::path dir=argv[1]; const auto seed=std::stoull(argv[2]);const int repeats=std::stoi(argv[3]);
    require(repeats>0 && repeats<=10,"invalid repetitions");
    std::filesystem::create_directories(dir);
    require(!std::filesystem::exists(dir/"initial.tsv"),"refusing to overwrite prior experiment");
    Verbosity=0; TimeLimit=std::numeric_limits<time_t>::max();
    StructureOptimization_SaveConfigurationInterval=1000000000;
    StructureOptimization_SaveNumberedConfiguration=false;
    std::mt19937_64 rng(seed);std::uniform_real_distribution<double> uniform(0,1);
    Configuration initial=emptyConfig();
    for(size_t i=0;i<N;++i){GeometryVector s(2);s.x[0]=uniform(rng);s.x[1]=uniform(rng);initial.Insert("A",s);}
    save(initial,dir/"initial.tsv");initial=load(dir/"initial.tsv");
    auto proto=potential(initial,dir/"modes.tsv");
    std::cout<<std::setprecision(17)<<"{\"seed\":"<<seed<<",\"N\":"<<N<<",\"M\":"<<M
      <<",\"chi\":"<<double(M)/(2*(N-1))<<",\"initial_validation\":";
    print(validate(initial,proto,dir/"initial_structure.tsv"));
    std::cout<<",\"initial_microbench\":";microbench(initial,proto);std::cout<<",\"pairs\":[";
    bool success=true;
    for(int rep=0;rep<repeats;++rep) {
        if(rep)std::cout<<',';
        bool fastFirst=(seed+rep)%2;
        auto first=optimize(dir,proto,fastFirst?Backend::FinufftType1:Backend::Direct,rep);
        auto second=optimize(dir,proto,fastFirst?Backend::Direct:Backend::FinufftType1,rep);
        const Run& d=fastFirst?second:first;const Run& f=fastFirst?first:second;
        success=success && d.validation.ground && f.validation.ground;
        std::cout<<"{\"repeat\":"<<rep<<",\"first\":\""<<(fastFirst?"finufft":"direct")<<"\",\"direct\":";print(d);
        std::cout<<",\"finufft\":";print(f);std::cout<<",\"positions\":";displacement(d.config,f.config);
        if(rep==0) {std::cout<<",\"direct_endpoint_microbench\":";microbench(d.config,proto);
            std::cout<<",\"finufft_endpoint_microbench\":";microbench(f.config,proto);}
        std::cout<<'}'<<std::flush;
        std::cerr<<"seed="<<seed<<" repeat="<<rep<<" E_direct="<<d.validation.energy<<" E_finufft="<<f.validation.energy
          <<" seconds="<<d.wall<<','<<f.wall<<" accepted="<<(d.validation.ground&&f.validation.ground)<<'\n';
    }
    std::cout<<"],\"all_ground_states\":"<<(success?"true":"false")<<"}\n";
    return success?0:2;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
