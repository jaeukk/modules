#ifndef PARAMETERINPUT_INCLUDED
#define PARAMETERINPUT_INCLUDED

/** \file ParameterInput.h
 *  \brief Small stream-based parameter-input helpers split out of etc.h: the
 *         type-erased value holder `cont` and the named-parameter reader
 *         `Parameters`, plus the variable-base odometer `n_naryNumber`.
 *         None currently have callers in the tree. */

#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <algorithm>   // std::transform in Parameters::Insert
#include <cctype>      // ::tolower in Parameters::Insert

/** \brief simulate binary and other number systems; unit in each digit can vary */
struct n_naryNumber {
protected:
	std::vector<unsigned int> digits;
	std::vector<unsigned int> upperbounds;
public:
	/** \brief generate 0 of numDigits; n-th digit: 0 ~ upperbounds[n]-1
	 *  \param numDigits    Number of digits.
	 *  \param upperbounds  Per-digit upper bounds (size must equal numDigits). */
	n_naryNumber(unsigned int numDigits, const std::vector<unsigned int> & upperbounds);

	/** \brief Construct a number of numDigits digits each sharing the same upper bound n.
	 *  \param numDigits  Number of digits.
	 *  \param n          Common upper bound for every digit. */
	n_naryNumber(unsigned int numDigits, unsigned int n);

	/** \brief Value of the N-th digit.
	 *  \param N  Digit index.
	 *  \return   The digit value, or 0 if N is out of range. */
	int GetNthDigit(size_t N) const;
	/** \brief increment in the number
	 *  \return  Reference to this number after incrementing. */
	n_naryNumber& operator++ (int);
};


/** Data structures to help inserting parameters. **/
#define STR(name) std::string(#name)

/**\brief A sloppy container to deal with pointers of four types (int, size_t, double, string).*/
struct cont {
	/** \brief A constructor for a int pointer.
	 *  \param _val  The int variable whose address is stored. */
	cont(int & _val) : cont() { _int_ptr = &_val; }				//< A constructor for a int pointer.
	/** \brief A constructor for a double pointer.
	 *  \param _val  The double variable whose address is stored. */
	cont(double & _val) : cont() { _double_ptr = &_val; }		//< A constructor for a double pointer.
	/** \brief A constructor for a size_t pointer.
	 *  \param _val  The size_t variable whose address is stored. */
	cont(size_t & _val) : cont() { _size_ptr = &_val; }			//< A constructor for a size_t pointer.
	/** \brief A constructor for a string pointer.
	 *  \param _val  The string variable whose address is stored. */
	cont(std::string & _val) : cont() { _string_ptr = &_val; }	//< A constructor for a string pointer.
																/** A copy constructor.
	 *  \param source  The cont to copy from. */
	cont(const cont & source) : _int_ptr(source._int_ptr), _double_ptr(source._double_ptr), _size_ptr(source._size_ptr), _string_ptr(source._string_ptr) {
		Fail = false;
	}
	/** \return this->Fail. True = operator >> failed. */
	bool IsFail() { return Fail; }
	/** \return the type of the non-null pointer. */
	std::string GetType();

	friend inline std::istream & operator >> (std::istream & is, cont & a);
	friend inline std::ostream & operator << (std::ostream & os, cont & a);
protected:
	//< pointers of four different types. Only one should be used at once. 
	int * _int_ptr;
	double * _double_ptr;
	size_t * _size_ptr;
	std::string * _string_ptr;

	bool Fail;	//< True = operator >> failed.
private:
	/** \brief A default constructor. */
	//A default constructor.
	cont() {
		_int_ptr = nullptr;	_double_ptr = nullptr;	_size_ptr = nullptr;	_string_ptr = nullptr; Fail = false;
	}
};
/** \brief Extraction operator for the cont struct.
 *  \param is  The input stream to read from.
 *  \param a   The cont whose pointed-to value receives the extracted value.
 *  \return    The input stream. */
std::istream & operator >> (std::istream & is, cont & a);
/** \brief Insertion operator for the cont struct.
 *  \param os  The output stream to write to.
 *  \param a   The cont whose pointed-to value is written.
 *  \return    The output stream. */
std::ostream & operator << (std::ostream & os, cont & a);


/** \brief Data structure to help input a set of parameters. */
struct Parameters {
	/** \brief A default constructor
	 *  \param Case_Sensitive  true = parameter names are matched case-sensitively. */
	//A default constructor
	Parameters(bool Case_Sensitive = false) { case_insensitive = Case_Sensitive; }

	/** A member function to add one item in the parameter list
	 * \param param_name	Name of the parameter.
	 * \param val			A variable into which an input parameter is transferred. 
							The type should be one of int, double, size_t, or std::string.	 */
	template<typename T> void Insert(const std::string & param_name, T & val) {
		std::string new_param_name(param_name);
		if (!this->case_insensitive) {
			// Not Case insensitive -> change the parameter name in lower cases.
			std::transform(new_param_name.begin(), new_param_name.end(), new_param_name.begin(), ::tolower);
		}
		ParamList.emplace_back(new_param_name, cont(val));
	}

	/** \brief A routine to input parameters.
	 *	The routine stops when one inputs "Start", "start", or "START".
	 *	\param ifile	An inputstream. */
	void InputRoutine(std::istream & ifile);

	/** \brief Display all parameter names. 
	 *	\param ofile	An outputstream. */
	void ShowAllParameterNames(std::ostream & ofile);

	/** \brief Display all parameter names and their values. 
	 *	\param ofile	An outputstream. */
	void ShowAllParameters(std::ostream & ofile);

protected:
	std::vector<std::pair<std::string, cont >>	ParamList;	//< ParamList[i].first = Parameter name, and ParamList[i].second = a pointer 
	size_t limit_wrong_input = 4;	//< When there are incorrect parameters more than "limit_wrong_input", InputRoutine(...) automatically ends.
	bool case_insensitive;			//< true = all parameter names are in lower case.
};

#endif
