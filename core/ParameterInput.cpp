/** \file ParameterInput.cpp
 *  \brief Implementations for the helpers declared in ParameterInput.h
 *         (n_naryNumber, cont, Parameters), split out of etc.cpp. */
#include "ParameterInput.h"
#include <algorithm>
#include <cctype>

n_naryNumber::n_naryNumber(unsigned int numDigits, const std::vector<unsigned int> & upperbounds) {
	if (numDigits == upperbounds.size()) {
		this->upperbounds = std::vector<unsigned int>(upperbounds);
		this->digits = std::vector<unsigned int>(numDigits, 0);
	}
	else {
		std::cerr << "Dimensions are mismatched\n";
	}
}

n_naryNumber::n_naryNumber(unsigned int numDigits, unsigned int n) : n_naryNumber::n_naryNumber(numDigits, std::vector<unsigned int>(numDigits, n)) {
}

int n_naryNumber::GetNthDigit(size_t N) const {
	if (N < digits.size())
		return this->digits[N];
	else {
		std::cerr << "exceed the digits\n";
		return 0;
	}
}

n_naryNumber& n_naryNumber::operator++ (int) {
	this->digits[0]++;
	for (size_t i = 0; i < this->digits.size(); i++) {
		if (this->digits[i] >= this->upperbounds[i]) {
			this->digits[i] = 0;
			if (i + 1 == this->digits.size()) {
				std::cout << "Overflow occurs!\n";
			}
			else
				this->digits[i + 1]++;
		}
		else
			break;
	}


	return *this;
}

// Function implementations of the cont class

std::string cont::GetType() {
	if (_int_ptr != nullptr)
		return "int";
	if (_double_ptr != nullptr)
		return "double";
	if (_size_ptr != nullptr)
		return "size";
	if (_string_ptr != nullptr)
		return "string";
	else {
		std::cerr << "Type is undefined\n";
		return "";
	}
}

/** \brief Extraction operator for the cont struct.
 *  \param is  The input stream to read from.
 *  \param a   The cont whose pointed-to value receives the extracted value.
 *  \return    The input stream. */
inline std::istream & operator >> (std::istream & is, cont & a) {
	if (a._int_ptr != nullptr) {
		is >> *(a._int_ptr);
	}
	else if (a._double_ptr != nullptr) {
		is >> *(a._double_ptr);
	}
	else if (a._size_ptr != nullptr) {
		is >> *(a._size_ptr);
	}
	else if (a._string_ptr != nullptr) {
		is >> *(a._string_ptr);
	}

	if (is.fail()) {
		is.clear();
		std::string dump;
		is >> dump;
		a.Fail = true;
		std::cout << "Fail: " << dump << " is an inappropriate parameter. \n";
	}
	return is;
}

/** \brief Insertion operator for the cont struct.
 *  \param os  The output stream to write to.
 *  \param a   The cont whose pointed-to value is written.
 *  \return    The output stream. */
inline std::ostream & operator << (std::ostream & os, cont & a) {
	if (a._int_ptr != nullptr) {
		os << *(a._int_ptr);
	}
	else if (a._double_ptr != nullptr) {
		os << *(a._double_ptr);
	}
	else if (a._size_ptr != nullptr) {
		os << *(a._size_ptr);
	}
	else if (a._string_ptr != nullptr) {
		os << *(a._string_ptr);
	}
	return os;
}


/* Member functions of the Parameters class. */

/** \brief A routine to input parameters.
 *	The routine stops when one inputs "Start", "start", or "START".
 *	\param ifile	An inputstream. */
void Parameters::InputRoutine(std::istream & ifile) {
	std::string tempstring;
	size_t count_wrong_input = 0;
	bool correct_input = false;
	for (;;) {
		ifile >> tempstring;
		if (!case_insensitive)
			std::transform(tempstring.begin(), tempstring.end(), tempstring.begin(), ::tolower);

		//End the routine.
		if (tempstring.compare("start") == 0 || tempstring.compare("Start") == 0 || tempstring.compare("START") == 0) {
			return;
		}

		for (auto i = this->ParamList.begin(); i != this->ParamList.end(); i++) {
			if (i->first.compare(tempstring) == 0) {
				ifile >> i->second;
				correct_input |= !(i->second.IsFail());
				break;
			}

			if (i == this->ParamList.end() - 1)
				std::cerr << "There is no parameter name that corresponds to " << tempstring << "\n";
		}

		if (!correct_input) {
			count_wrong_input++;

			if (count_wrong_input >= limit_wrong_input) {
				std::cerr << " There were too many incorrect inputs!\n";
				return;
			}
		}

		correct_input = false;
	}
}

/** \brief Display all parameter names.
 *	\param ofile	An outputstream. */
void Parameters::ShowAllParameterNames(std::ostream & ofile) {
	ofile << "------------------------------------------\n";
	ofile << "\t\tList of parameter names\n";
	ofile << "------------------------------------------\n";
	for (auto i = ParamList.begin(); i != ParamList.end(); i++) {
		ofile << i->first << "(" << i->second.GetType() << ")\n";
	}
}

/** \brief Display all parameter names and their values.
 *	\param ofile	An outputstream. */
void Parameters::ShowAllParameters(std::ostream & ofile) {
	ofile << "------------------------------------------\n";
	ofile << "\t\tList of parameters\n";
	ofile << "------------------------------------------\n";
	for (auto i = ParamList.begin(); i != ParamList.end(); i++)
		ofile << i->first << " :\t\t" << i->second << "\n";
}
