/**************************************************************************

   Copyright © 2018 Sergey Poluyan <svpoluyan@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

**************************************************************************/
#include <numeric/NumericTraits.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <core/conformation/util.hh>

#include "bbtools.hh"

namespace pepsgo
{
namespace bbtools
{
    
core::Real to_positive_radians(core::Real radian)
{
//    if(radian > 0)
//        return radian;
//    numeric::Real delta = numeric::NumericTraits< core::Real >::pi() - std::abs(radian);
//    return (numeric::NumericTraits< core::Real >::pi() + delta);
	if (radian > 0)
		return radian;
	while (radian < 0)
		radian += 2.0 * numeric::NumericTraits < core::Real > ::pi();
	return radian;
}
core::Real to_negative_radians(core::Real radian)
{
	if (radian < 0)
		return radian;
	while (radian > 0)
		radian -= 2.0 * numeric::NumericTraits<core::Real>::pi();
	return radian;
}    
core::Real normalize_to_mpi_to_ppi(core::Real radian)
{
	if (radian < numeric::NumericTraits<core::Real>::pi() && radian > -numeric::NumericTraits<core::Real>::pi())
		return radian;
	if (radian > numeric::NumericTraits<core::Real>::pi())
		return radian - 2.0 * numeric::NumericTraits<core::Real>::pi();
	else
		return radian + 2.0 * numeric::NumericTraits<core::Real>::pi();
}

void make_ideal_peptide(core::pose::Pose& ideal_peptide, const core::pose::Pose& peptide)
{
	ideal_peptide = peptide;
	core::conformation::Conformation conf = ideal_peptide.conformation();

	for (numeric::Size seqpos = 1, seqpos_end = ideal_peptide.total_residue(); seqpos <= seqpos_end; ++seqpos)
	{
		core::conformation::idealize_position(seqpos, ideal_peptide.conformation());
	}
}

}
}