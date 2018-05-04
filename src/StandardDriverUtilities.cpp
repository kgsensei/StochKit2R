/*!

*/
#include "StandardDriverUtilities.h"

namespace STOCHKIT
{


	std::string StandardDriverUtilities::size_t2string(std::size_t number) {
		std::string theString;
		std::stringstream ss;
		ss<<number;
		theString=ss.str();
		return theString;
	}

}
