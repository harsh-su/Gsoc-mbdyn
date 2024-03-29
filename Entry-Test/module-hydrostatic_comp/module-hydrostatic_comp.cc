/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "constltp.h"

class HydrostaticCompressionCL
: public ConstitutiveLaw<Vec3,Mat3x3> {
private:
	doublereal m_K;
	doublereal hydrostatic_pressure;
#ifdef USE_NETCDF
	MBDynNcVar Var_HydrostaticP;
#endif // USE_NETCDF
public:
	HydrostaticCompressionCL(doublereal K)
	: m_K(K){
		ConstitutiveLaw<Vec3, Mat3x3>::FDE= mb_deye<Mat3x3>(m_K);
	};

	virtual ~HydrostaticCompressionCL(void) {
		NO_OP;
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual ConstitutiveLaw<Vec3,Mat3x3>* pCopy(void) const {
		ConstitutiveLaw<Vec3,Mat3x3>* pCL = 0;

		typedef HydrostaticCompressionCL cl;
		SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(m_K));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "Modulus_of_compressiblity, " << m_K ;
	};

	virtual void Update(const Vec3& Eps , const Vec3& Epsilon) {
		ConstitutiveLaw<Vec3,Mat3x3>::Epsilon = Eps;
        ConstitutiveLaw<Vec3,Mat3x3>::F = ConstitutiveLaw<Vec3,Mat3x3>::Epsilon*m_K;
		/* Now to Calculate hydrostatic Pressure */
		hydrostatic_pressure = -(F[0,0] + F[1,1] + F[2,2]) / 3.;
        
	};

	virtual std::ostream& OutputAppend(std::ostream& out) const{
		return out << " " << hydrostatic_pressure;
	};


	virtual void NetCDFOutputAppend(OutputHandler& OH) const
	{
	#ifdef USE_NETCDF
	OH.WriteNcVar(Var_HydrostaticP, hydrostatic_pressure);
	#endif // USE_NETCDF
	};



	virtual void OutputAppendPrepare(OutputHandler& OH, const std::string& name){
	#ifdef USE_NETCDF
	ASSERT(OH.IsOpen(OutputHandler::NETCDF));
	if (OH.UseNetCDF(OutputHandler::LOADABLE)) 
	{
		Var_HydrostaticP = OH.CreateVar<doublereal>(name + ".a", 
				OutputHandler::Dimensions::Dimensionless, 
				"Hydrostatic_Pressure");
		
	}
	#endif // USE_NETCDF
	};	
};

/* specific functional object(s) */
struct HydrostaticCompressionCLR : public ConstitutiveLawRead<Vec3,Mat3x3> {
	virtual ConstitutiveLaw<Vec3,Mat3x3> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<Vec3,Mat3x3>* pCL = 0;

		CLType = ConstLawType::ELASTIC;

		doublereal dS = HP.GetReal();
		if (dS == 0.) {
			silent_cerr("warning, null value " << HP.GetLineData() << std::endl);
		}

		typedef HydrostaticCompressionCL L;
		SAFENEWWITHCONSTRUCTOR(pCL, L, L(dS));

		return pCL;
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	ConstitutiveLawRead<Vec3, Mat3x3> *rf3D = new HydrostaticCompressionCLR;
	if (!SetCL3D("HydrostaticCompression", rf3D)) {
		delete rf3D;

		silent_cerr("HydrostaticCompressionConstitutiveLaw3D: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}