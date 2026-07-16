#include "Trajectory.h"
#include <sstream>

Trajectory::Trajectory(const Trajectory & src){
	this->dimension = src.dimension;
	this->sigma_host = src.sigma_host; this->sigma_particle = src.sigma_particle;
	this->R2_host = src.R2_host;	this->R2_particle = src.R2_particle;
	this->R2_s = src.R2_s; tau_s = src.tau_s;
	this->tau_total = src.tau_total;

	this->pos_curr = GeometryVector(src.pos_curr);
	this->delta_pos = GeometryVector(src.delta_pos);
}

void Trajectory::ReadLog(std::istream & ifile){
	std::string phrase("#d");
	size_t dim = 0;

	std::string buffer; std::istringstream iss;

	while (dim == 0 && std::getline(ifile, buffer)){
		if(buffer[0] == '#'){
			dim = (size_t) std::stoi(buffer.substr(phrase.size()));	
		}
	}

	if (dim == this->dimension){
		// Read pos_curr
		std::getline(ifile, buffer);
		this->pos_curr.SetDimension(dim);
		iss.str(buffer);
		iss.clear();
		for (size_t j=0; j<dim; j++)
			iss >> this->pos_curr.x[j];
		
//		if (ifile.eof() == false)
//			vertices.push_back(temp);

		// Read delta_pos
		std::getline(ifile, buffer);
		this->delta_pos.SetDimension(dim);
		iss.str(buffer);
		iss.clear();
		for (size_t j=0; j<dim; j++)
			iss >> this->delta_pos.x[j];


		std::getline(ifile, buffer);
		iss.str(buffer);
		iss.clear();
		double temp1 = 0.0, temp2 = 0.0;
		// Read R2_host
		iss >> temp1;
		this->Update_R2_host(temp1);

		// Read R2_particle
		iss >> temp1;
		this->Update_R2_particle(temp1);

		// Read R2_s, tau_s
		iss >> temp1; iss >> temp2;
		this->Update_R2_s(temp1, temp2);
		
		// compare time
		iss >> temp1;
		double err_time = temp1 - this->GetTimeTotal();
		this->tau_total += err_time;
	}
	else{
		std::cerr<<"Trajectory::ReadLog(...)\twrong dimension!\n";
	}
}

std::ostream & operator << (std::ostream & os, const Trajectory & t) 
{
//	os << "[Trajectory]\n";
	os << "#d\t" << t.dimension<<"\n";
	os << t.pos_curr;
	os << t.delta_pos;

	os << t.R2_host << "\t" << t.R2_particle <<"\t";
	os << t.R2_s << "\t" << t.tau_s << "\t";
	os << t.tau_total<<"\n";
	return os;
}
