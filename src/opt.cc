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

#include <opt.hh>

namespace pepsgo
{

opt_element::opt_element()
{

}

opt_element::opt_element(core::id::DOF_ID _dofid,
                         std::pair<core::Real, core::Real> _bounds,
                         core::chemical::AA _amino_acid,
                         std::string _torsion_name,
                         core::Size _chainid,
                         core::Size _seqpos,
                         core::Size _nchi,
                         core::Size _ichi)
{
    dofid = _dofid;
    bounds = _bounds;
    amino_acid = _amino_acid;
    torsion_name = _torsion_name;
    chainid = _chainid;
    seqpos = _seqpos;
    nchi = _nchi;
    ichi = _ichi;
}


}
