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
#ifndef INCLUDED_opt_hh
#define INCLUDED_opt_hh

#include <core/pose/Pose.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/TorsionID.hh>

#include <tuple>

namespace pepsgo
{

struct opt_element
{
    core::id::DOF_ID dofid;
    std::pair<core::Real, core::Real> bounds;
    core::chemical::AA amino_acid;
    std::string torsion_name;
    core::Size chainid;
    core::Size seqpos;
    core::Size nchi;
    core::Size ichi;

    opt_element();

    opt_element(core::id::DOF_ID _dofid, std::pair<core::Real, core::Real> _bounds, core::chemical::AA _amino_acid, std::string _torsion_name, core::Size _chainid,
                core::Size _seqpos,
                core::Size _nchi,
                core::Size _ichi);

    friend std::ostream& operator<<(std::ostream& stream, const opt_element& obj)
    {
        stream << obj.chainid << ' ' << obj.amino_acid << ' ' << obj.torsion_name << ' ' << obj.seqpos << ' ' << obj.nchi << ' ' << obj.ichi << ' ' << obj.bounds.first << ' ' << obj.bounds.second;
        return stream;
    }
};

struct ranges
{
    std::tuple<bool, size_t, size_t> phipsi;
    std::tuple<bool, size_t, size_t> omega;
    std::tuple<bool, size_t, size_t> chi;

    bool do_phipsi;
    bool do_omega;
    bool do_chi;
};

}


#endif
