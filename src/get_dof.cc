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
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/util.hh>

#include <numeric/NumericTraits.hh>

#include <get_dof.hh>

namespace pepsgo
{

std::vector<core::id::DOF_ID> get_phi_psi(core::pose::Pose& pose)
{
    std::vector < core::id::DOF_ID > peptide_dof_BB;
    for(core::Size seqpos = 1, seqpos_end = pose.total_residue(); seqpos <= seqpos_end; ++seqpos)
    {
        core::id::TorsionID bb_phi_tid(seqpos, core::id::BB, core::id::phi_torsion);
        core::id::TorsionID bb_psi_tid(seqpos, core::id::BB, core::id::psi_torsion);

        peptide_dof_BB.push_back(pose.conformation().dof_id_from_torsion_id(bb_phi_tid));
        peptide_dof_BB.push_back(pose.conformation().dof_id_from_torsion_id(bb_psi_tid));
    }

    peptide_dof_BB.erase(peptide_dof_BB.begin());
    peptide_dof_BB.erase(peptide_dof_BB.end());

    return peptide_dof_BB;
}

std::vector<opt_element> get_phi_psi_with_info(core::pose::Pose& pose)
{
    std::vector<opt_element> peptide_dof_BB;
    for(core::Size seqpos = 1, seqpos_end = pose.total_residue(); seqpos <= seqpos_end; ++seqpos)
    {
        core::id::TorsionID bb_phi_tid(seqpos, core::id::BB, core::id::phi_torsion);
        core::id::TorsionID bb_psi_tid(seqpos, core::id::BB, core::id::psi_torsion);

        peptide_dof_BB.push_back(
            opt_element(pose.conformation().dof_id_from_torsion_id(bb_phi_tid), std::make_pair(
                            -numeric::NumericTraits<core::Real>::pi(), numeric::NumericTraits<core::Real>::pi()
                        ), pose.residue(seqpos).type().aa(), "bbphi", 0, seqpos,
                        0, 0));
        peptide_dof_BB.push_back(
            opt_element(pose.conformation().dof_id_from_torsion_id(bb_psi_tid),
                        std::make_pair(-numeric::NumericTraits<core::Real>::pi(), numeric::NumericTraits<core::Real>::pi()), pose.residue(seqpos).type().aa(), "bbpsi", 0, seqpos,
                        0,
                        0));
    }

    peptide_dof_BB.erase(peptide_dof_BB.begin());
    peptide_dof_BB.erase(peptide_dof_BB.end());

    return peptide_dof_BB;
}
std::vector<opt_element> get_phi_psi_with_info(core::pose::Pose& pose, std::vector<core::Size> &ind, size_t chain_id, bool remove_last)
{
    std::vector<opt_element> peptide_dof_BB;
    for(core::Size seqpos = 0, seqpos_end = ind.size(); seqpos != seqpos_end; ++seqpos)
    {
        core::id::TorsionID bb_phi_tid(ind[seqpos], core::id::BB, core::id::phi_torsion);
        core::id::TorsionID bb_psi_tid(ind[seqpos], core::id::BB, core::id::psi_torsion);

        peptide_dof_BB.push_back(
            opt_element(pose.conformation().dof_id_from_torsion_id(bb_phi_tid),
                        std::make_pair(-numeric::NumericTraits<core::Real>::pi(), numeric::NumericTraits<core::Real>::pi()), pose.residue(ind[seqpos]).type().aa(), "bbphi",
                        chain_id, ind[seqpos], 0, 0));
        peptide_dof_BB.push_back(
            opt_element(pose.conformation().dof_id_from_torsion_id(bb_psi_tid),
                        std::make_pair(-numeric::NumericTraits<core::Real>::pi(), numeric::NumericTraits<core::Real>::pi()), pose.residue(ind[seqpos]).type().aa(), "bbpsi",
                        chain_id, ind[seqpos], 0, 0));
    }

    if(remove_last)
    {
        peptide_dof_BB.erase(peptide_dof_BB.begin());
        peptide_dof_BB.erase(peptide_dof_BB.end());
    }
    return peptide_dof_BB;
}


std::vector<core::id::DOF_ID> get_peptide_sc_dof_THETA(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector < core::id::DOF_ID > peptide_sc_dof_THETA;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 1, atom_end = rt->all_sc_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->all_sc_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(atom_type == 1)
            {
                if(rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                {
                    peptide_sc_dof_THETA.push_back(core::id::DOF_ID(aid, core::id::THETA));
                }
            }
            if(atom_type == 2)
            {
                if(!rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                {
                    if(rt->atom_type(rt->all_sc_atoms().at(atom)).is_aromatic())
                    {
                        continue;
                    }
                    peptide_sc_dof_THETA.push_back(core::id::DOF_ID(aid, core::id::THETA));
                }
            }
            if(atom_type == 0)
            {
                peptide_sc_dof_THETA.push_back(core::id::DOF_ID(aid, core::id::THETA));
            }
        }
    }
    return peptide_sc_dof_THETA;
}

std::vector<opt_element> get_peptide_sc_dof_THETA_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector<opt_element> peptide_sc_dof_THETA;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 1, atom_end = rt->all_sc_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->all_sc_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(atom_type == 1)
            {
                if(rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                {
                    peptide_sc_dof_THETA.push_back(

                        opt_element(core::id::DOF_ID(aid, core::id::THETA), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "sc theta hydrogen", 0, seqpos, 0, 0)

                    );
                }
            }
            if(atom_type == 2)
            {
                if(!rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                {
                    if(rt->atom_type(rt->all_sc_atoms().at(atom)).is_aromatic())
                    {
                        continue;
                    }
                    peptide_sc_dof_THETA.push_back(

                        opt_element(core::id::DOF_ID(aid, core::id::THETA), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "sc theta nothydrogen", 0, seqpos, 0, 0)

                    );
                }
            }
            if(atom_type == 0)
            {
                peptide_sc_dof_THETA.push_back(

                    opt_element(core::id::DOF_ID(aid, core::id::THETA), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "sc theta", 0, seqpos, 0, 0)

                );
            }
        }
    }
    return peptide_sc_dof_THETA;
}


std::vector<core::id::DOF_ID> get_peptide_sc_dof_D(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector < core::id::DOF_ID > peptide_sc_dof_D;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 1, atom_end = rt->all_sc_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->all_sc_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(atom_type == 1)
            {
                if(rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                {
                    peptide_sc_dof_D.push_back(core::id::DOF_ID(aid, core::id::D));
                }
            }
            if(atom_type == 2)
            {
                if(!rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                {
                    if(rt->atom_type(rt->all_sc_atoms().at(atom)).is_aromatic())
                    {
                        continue;
                    }
                    peptide_sc_dof_D.push_back(core::id::DOF_ID(aid, core::id::D));
                }
            }
            if(atom_type == 0)
            {
                peptide_sc_dof_D.push_back(core::id::DOF_ID(aid, core::id::D));
            }
        }
    }
    return peptide_sc_dof_D;
}

std::vector<opt_element> get_peptide_sc_dof_D_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector<opt_element> peptide_sc_dof_D;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 1, atom_end = rt->all_sc_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->all_sc_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(atom_type == 1)
            {
                if(rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                {
                    peptide_sc_dof_D.push_back(

                        opt_element(core::id::DOF_ID(aid, core::id::D), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "sc d hydrogen", 0, seqpos, 0, 0)

                    );
                }
            }
            if(atom_type == 2)
            {
                if(!rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                {
                    if(rt->atom_type(rt->all_sc_atoms().at(atom)).is_aromatic())
                    {
                        continue;
                    }
                    peptide_sc_dof_D.push_back(
                        opt_element(core::id::DOF_ID(aid, core::id::D), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "sc d nothydrogen", 0, seqpos, 0, 0));
                }
            }
            if(atom_type == 0)
            {
                peptide_sc_dof_D.push_back(

                    opt_element(core::id::DOF_ID(aid, core::id::D), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "sc d", 0, seqpos, 0, 0)

                );
            }
        }
    }
    return peptide_sc_dof_D;
}

std::vector<core::id::DOF_ID> get_peptide_sc_dof_PHI_without_CHI(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector < core::id::DOF_ID > peptide_sc_dof_PHI;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));

        utility::vector1 < core::id::DOF_ID > chi_tids(pose.residue(seqpos).nchi());
        if(pose.residue(seqpos).nchi())
        {
            for(core::Size i = 1, n = pose.residue(seqpos).nchi(); i <= n; ++i)
                chi_tids[i] = pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, i));
        }

        for(core::Size atom = 1, atom_end = rt->all_sc_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            bool flag_add = true;
            core::id::AtomID aid(rt->all_sc_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            core::id::DOF_ID dof_PHI(aid, core::id::PHI);
            if(pose.residue(seqpos).nchi())
            {
                for(core::Size i = 1, n = pose.residue(seqpos).nchi(); i <= n; ++i)
                {
                    if(chi_tids[i] == dof_PHI)
                    {
                        flag_add = false;
                        break;
                    }
                }
            }
            if(flag_add)
            {
                if(atom_type == 1)
                {
                    if(rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                    {
                        peptide_sc_dof_PHI.push_back(core::id::DOF_ID(aid, core::id::PHI));
                    }
                }
                if(atom_type == 2)
                {
                    if(!rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                    {

                        /*if(std::abs(
                         std::abs( pose.dof(core::id::DOF_ID(aid, core::id::PHI)) ) -
                         numeric::NumericTraits< core::Real >::pi() )
                         < 0.01 )
                         continue;
                         if(std::abs(
                         std::abs( pose.dof(core::id::DOF_ID(aid, core::id::PHI)) ) )
                         < 0.01 )
                         continue;*/

                        if(rt->atom_type(rt->all_sc_atoms().at(atom)).is_aromatic())
                        {
                            continue;
                        }

                        peptide_sc_dof_PHI.push_back(core::id::DOF_ID(aid, core::id::PHI));
                    }
                }
                if(atom_type == 0)
                {
                    peptide_sc_dof_PHI.push_back(core::id::DOF_ID(aid, core::id::PHI));
                }
            }
        }
    }
    return peptide_sc_dof_PHI;
}
std::vector<opt_element> get_peptide_sc_dof_PHI_without_CHI_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector<opt_element> peptide_sc_dof_PHI;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));

        utility::vector1<core::id::DOF_ID> chi_tids(pose.residue(seqpos).nchi());
        if(pose.residue(seqpos).nchi())
        {
            for(core::Size i = 1, n = pose.residue(seqpos).nchi(); i <= n; ++i)
                chi_tids[i] = pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, i));
        }

        for(core::Size atom = 1, atom_end = rt->all_sc_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            bool flag_add = true;
            core::id::AtomID aid(rt->all_sc_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            core::id::DOF_ID dof_PHI(aid, core::id::PHI);
            if(pose.residue(seqpos).nchi())
            {
                for(core::Size i = 1, n = pose.residue(seqpos).nchi(); i <= n; ++i)
                {
                    if(chi_tids[i] == dof_PHI)
                    {
                        flag_add = false;
                        break;
                    }
                }
            }
            if(flag_add)
            {
                if(atom_type == 1)
                {
                    if(rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                    {
                        peptide_sc_dof_PHI.push_back(
                            opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "sc phi hydrogen", 0, seqpos, 0, 0));
                    }
                }
                if(atom_type == 2)
                {
                    if(!rt->atom(rt->all_sc_atoms().at(atom)).is_hydrogen())
                    {

                        /*if(std::abs(
                         std::abs( pose.dof(core::id::DOF_ID(aid, core::id::PHI)) ) -
                         numeric::NumericTraits< core::Real >::pi() )
                         < 0.01 )
                         continue;
                         if(std::abs(
                         std::abs( pose.dof(core::id::DOF_ID(aid, core::id::PHI)) ) )
                         < 0.01 )
                         continue;*/

                        if(rt->atom_type(rt->all_sc_atoms().at(atom)).is_aromatic())
                        {
                            continue;
                        }

                        peptide_sc_dof_PHI.push_back(

                            opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "sc phi nothydrogen", 0, seqpos, 0, 0)

                        );
                    }
                }
                if(atom_type == 0)
                {
                    peptide_sc_dof_PHI.push_back(

                        opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "sc phi", 0, seqpos, 0, 0)

                    );
                }
            }
        }
    }
    return peptide_sc_dof_PHI;
}
std::vector<core::id::DOF_ID> get_peptide_bb_dof_THETA_without_mc(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector < core::id::DOF_ID > peptide_bb_dof_THETA_without_mc;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 4, atom_end = rt->all_bb_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->all_bb_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(atom_type == 1)
            {
                if(rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_THETA_without_mc.push_back(core::id::DOF_ID(aid, core::id::THETA));
                }
            }
            if(atom_type == 2)
            {
                if(!rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_THETA_without_mc.push_back(core::id::DOF_ID(aid, core::id::THETA));
                }
            }
            if(atom_type == 0)
            {
                peptide_bb_dof_THETA_without_mc.push_back(core::id::DOF_ID(aid, core::id::THETA));
            }
        }
    }
    return peptide_bb_dof_THETA_without_mc;
}

std::vector<opt_element> get_peptide_bb_dof_THETA_without_mc_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector<opt_element> peptide_bb_dof_THETA_without_mc;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 4, atom_end = rt->all_bb_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->all_bb_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(atom_type == 1)
            {
                if(rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_THETA_without_mc.push_back(

                        opt_element(core::id::DOF_ID(aid, core::id::THETA), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "bb theta hydrogen", 0, seqpos, 0, 0)

                    );


                }
            }
            if(atom_type == 2)
            {
                if(!rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_THETA_without_mc.push_back(

                        opt_element(core::id::DOF_ID(aid, core::id::THETA), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "bb theta nothydrogen", 0, seqpos, 0, 0)

                    );
                }
            }
            if(atom_type == 0)
            {
                peptide_bb_dof_THETA_without_mc.push_back(
                    opt_element(core::id::DOF_ID(aid, core::id::THETA), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "bb theta", 0, seqpos, 0, 0)

                );
            }
        }
    }
    return peptide_bb_dof_THETA_without_mc;
}

std::vector<core::id::DOF_ID> get_peptide_bb_dof_D_without_mc(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector < core::id::DOF_ID > peptide_bb_dof_D_without_mc;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 4, atom_end = rt->all_bb_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->all_bb_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(atom_type == 1)
            {
                if(rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_D_without_mc.push_back(core::id::DOF_ID(aid, core::id::D));
                }
            }
            if(atom_type == 2)
            {
                if(!rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_D_without_mc.push_back(core::id::DOF_ID(aid, core::id::D));
                }
            }
            if(atom_type == 0)
            {
                peptide_bb_dof_D_without_mc.push_back(core::id::DOF_ID(aid, core::id::D));
            }
        }
    }
    return peptide_bb_dof_D_without_mc;
}

std::vector<opt_element> get_peptide_bb_dof_D_without_mc_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector<opt_element> peptide_bb_dof_D_without_mc;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 4, atom_end = rt->all_bb_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->all_bb_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(atom_type == 1)
            {
                if(rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_D_without_mc.push_back(

                        opt_element(core::id::DOF_ID(aid, core::id::D), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "bb d hydrogen", 0, seqpos, 0, 0));
                }
            }
            if(atom_type == 2)
            {
                if(!rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_D_without_mc.push_back(

                        opt_element(core::id::DOF_ID(aid, core::id::D), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "bb d nothydrogen", 0, seqpos, 0, 0)

                    );
                }
            }
            if(atom_type == 0)
            {
                peptide_bb_dof_D_without_mc.push_back(

                    opt_element(core::id::DOF_ID(aid, core::id::D), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "bb d nothydrogen", 0, seqpos, 0, 0)

                );
            }
        }
    }
    return peptide_bb_dof_D_without_mc;
}

std::vector<core::id::DOF_ID> get_peptide_bb_dof_PHI_without_mc_without_term(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector < core::id::DOF_ID > peptide_bb_dof_PHI_without_mc;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        //if(pose.residue(seqpos).is_terminus()) pose.residue(seqpos).is_lower_terminus() pose.residue(seqpos).is_upper_terminus()
        if(seqpos == 1 || seqpos == pose.total_residue())
        {
            continue;
        }
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 4, atom_end = rt->all_bb_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->all_bb_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            //std::cout << pose.residue(seqpos).atom_name(rt->all_bb_atoms().at(atom)) << std::endl;
            if(atom_type == 1)
            {
                if(rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_PHI_without_mc.push_back(core::id::DOF_ID(aid, core::id::PHI));
                }
            }
            if(atom_type == 2)
            {
                if(!rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_PHI_without_mc.push_back(core::id::DOF_ID(aid, core::id::PHI));
                }
            }
            if(atom_type == 0)
            {
                peptide_bb_dof_PHI_without_mc.push_back(core::id::DOF_ID(aid, core::id::PHI));
            }
        }
    }
    return peptide_bb_dof_PHI_without_mc;
}

std::vector<opt_element> get_peptide_bb_dof_PHI_without_mc_without_term_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type)
{
    std::vector<opt_element> peptide_bb_dof_PHI_without_mc;
    for(core::Size seqpos = start_ind, seqpos_end = stop_ind; seqpos <= seqpos_end; ++seqpos)
    {
        //if(pose.residue(seqpos).is_terminus()) pose.residue(seqpos).is_lower_terminus() pose.residue(seqpos).is_upper_terminus()
        if(seqpos == 1 || seqpos == pose.total_residue())
        {
            continue;
        }
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 4, atom_end = rt->all_bb_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->all_bb_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            //std::cout << pose.residue(seqpos).atom_name(rt->all_bb_atoms().at(atom)) << std::endl;
            if(atom_type == 1)
            {
                if(rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_PHI_without_mc.push_back(

                        opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "bb phi hydrogen", 0, seqpos, 0, 0)

                    );
                }
            }
            if(atom_type == 2)
            {
                if(!rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
                {
                    peptide_bb_dof_PHI_without_mc.push_back(
                        opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "bb phi nothydrogen", 0, seqpos, 0, 0)

                    );
                }
            }
            if(atom_type == 0)
            {
                peptide_bb_dof_PHI_without_mc.push_back(

                    opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "bb phi hydrogen", 0, seqpos, 0, 0)

                );
            }
        }
    }
    return peptide_bb_dof_PHI_without_mc;
}

std::vector<core::id::DOF_ID> get_peptide_bb_dof_PHI_first_term(core::pose::Pose& pose, int atom_type)
{
    std::vector < core::id::DOF_ID > peptide_bb_dof_PHI_term;
    core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(1).type()));
    for(core::Size atom = 4, atom_end = rt->all_bb_atoms().vector().size(); atom <= atom_end; ++atom)
    {
        if(atom == 5)
            continue;
        core::id::AtomID aid(rt->all_bb_atoms().at(atom), 1);
        if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
        {
            continue;
        }
        //std::cout << pose.residue(1).atom_name(rt->all_bb_atoms().at(atom)) << std::endl;
        if(atom_type == 1)
        {
            if(rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
            {
                peptide_bb_dof_PHI_term.push_back(core::id::DOF_ID(aid, core::id::PHI));
            }
        }
        if(atom_type == 2)
        {
            if(!rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
            {
                peptide_bb_dof_PHI_term.push_back(core::id::DOF_ID(aid, core::id::PHI));
            }
        }
        if(atom_type == 0)
        {
            peptide_bb_dof_PHI_term.push_back(core::id::DOF_ID(aid, core::id::PHI));
        }
    }
    return peptide_bb_dof_PHI_term;
}
std::vector<opt_element> get_peptide_bb_dof_PHI_first_term_with_info(core::pose::Pose& pose, int atom_type)
{
    std::vector<opt_element> peptide_bb_dof_PHI_term;
    core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(1).type()));
    for(core::Size atom = 4, atom_end = rt->all_bb_atoms().vector().size(); atom <= atom_end; ++atom)
    {
        if(atom == 5)
            continue;
        core::id::AtomID aid(rt->all_bb_atoms().at(atom), 1);
        if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
        {
            continue;
        }
        //std::cout << pose.residue(1).atom_name(rt->all_bb_atoms().at(atom)) << std::endl;
        if(atom_type == 1)
        {
            if(rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
            {
                peptide_bb_dof_PHI_term.push_back(

                    opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(1).type().aa(), "bb phi hydrogen", 0, 1, 0, 0)

                );
            }
        }
        if(atom_type == 2)
        {
            if(!rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
            {
                peptide_bb_dof_PHI_term.push_back(

                    opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(1).type().aa(), "bb phi hydrogen", 0, 1, 0, 0)

                );
            }
        }
        if(atom_type == 0)
        {
            peptide_bb_dof_PHI_term.push_back(

                opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(1).type().aa(), "bb phi hydrogen", 0, 1, 0, 0)

            );
        }
    }
    return peptide_bb_dof_PHI_term;
}

std::vector<core::id::DOF_ID> get_peptide_bb_dof_PHI_last_term(core::pose::Pose& pose, int atom_type)
{
    std::vector < core::id::DOF_ID > peptide_bb_dof_PHI_term;
    core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(pose.total_residue()).type()));
    for(core::Size atom = 5, atom_end = rt->all_bb_atoms().vector().size(); atom <= atom_end; ++atom)
    {
        core::id::AtomID aid(rt->all_bb_atoms().at(atom), pose.total_residue());
        if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
        {
            continue;
        }
        //std::cout << pose.residue(pose.total_residue()).atom_name(rt->all_bb_atoms().at(atom)) << std::endl;
        if(atom_type == 1)
        {
            if(rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
            {
                peptide_bb_dof_PHI_term.push_back(core::id::DOF_ID(aid, core::id::PHI));
            }
        }
        if(atom_type == 2)
        {
            if(!rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
            {
                peptide_bb_dof_PHI_term.push_back(core::id::DOF_ID(aid, core::id::PHI));
            }
        }
        if(atom_type == 0)
        {
            peptide_bb_dof_PHI_term.push_back(core::id::DOF_ID(aid, core::id::PHI));
        }
    }
    return peptide_bb_dof_PHI_term;
}

std::vector<opt_element> get_peptide_bb_dof_PHI_last_term_with_info(core::pose::Pose& pose, int atom_type)
{
    std::vector<opt_element> peptide_bb_dof_PHI_term;
    core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(pose.total_residue()).type()));
    for(core::Size atom = 5, atom_end = rt->all_bb_atoms().vector().size(); atom <= atom_end; ++atom)
    {
        core::id::AtomID aid(rt->all_bb_atoms().at(atom), pose.total_residue());
        if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
        {
            continue;
        }
        //std::cout << pose.residue(pose.total_residue()).atom_name(rt->all_bb_atoms().at(atom)) << std::endl;
        if(atom_type == 1)
        {
            if(rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
            {
                peptide_bb_dof_PHI_term.push_back(

                    opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(pose.total_residue()).type().aa(), "bb phi hydrogen",
                                0,
                                pose.total_residue(), 0, 0)

                );
            }
        }
        if(atom_type == 2)
        {
            if(!rt->atom(rt->all_bb_atoms().at(atom)).is_hydrogen())
            {
                peptide_bb_dof_PHI_term.push_back(

                    opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(pose.total_residue()).type().aa(), "bb phi hydrogen",
                                0,
                                pose.total_residue(), 0, 0)

                );
            }
        }
        if(atom_type == 0)
        {
            peptide_bb_dof_PHI_term.push_back(

                opt_element(core::id::DOF_ID(aid, core::id::PHI), std::make_pair(0, 0), pose.residue(pose.total_residue()).type().aa(), "bb phi hydrogen", 0,
                            pose.total_residue(), 0, 0)

            );
        }
    }
    return peptide_bb_dof_PHI_term;
}


std::vector<core::id::DOF_ID> get_peptide_bb_dof_PHI_term_chi1(core::pose::Pose& pose)
{
    std::vector < core::id::DOF_ID > peptide_bb_dof_PHI_term;
    core::chemical::ResidueTypeOP rt1(new core::chemical::ResidueType(pose.residue(1).type()));
    core::id::AtomID aid1(rt1->all_bb_atoms().at(5), 1);
    std::cout << pose.residue(1).atom_name(pose.residue(1).type().all_bb_atoms().at(5)) << std::endl;
    peptide_bb_dof_PHI_term.push_back(core::id::DOF_ID(aid1, core::id::PHI));
    return peptide_bb_dof_PHI_term;
}
std::vector<opt_element> get_peptide_bb_dof_PHI_term_chi1_with_info(core::pose::Pose& pose)
{
    std::vector<opt_element> peptide_bb_dof_PHI_term;
    core::chemical::ResidueTypeOP rt1(new core::chemical::ResidueType(pose.residue(1).type()));
    core::id::AtomID aid1(rt1->all_bb_atoms().at(5), 1);
    std::cout << pose.residue(1).atom_name(pose.residue(1).type().all_bb_atoms().at(5)) << std::endl;
    peptide_bb_dof_PHI_term.push_back(

        opt_element(core::id::DOF_ID(aid1, core::id::PHI), std::make_pair(0, 0), pose.residue(1).type().aa(), "bb phi chi1", 0, 1, 0, 0)

    );
    return peptide_bb_dof_PHI_term;
}

std::vector<core::id::DOF_ID> get_peptide_bb_dof_PHI_term_chi2(core::pose::Pose& pose)
{
    std::vector < core::id::DOF_ID > peptide_bb_dof_PHI_term;
    core::chemical::ResidueTypeOP rt2(new core::chemical::ResidueType(pose.residue(pose.total_residue()).type()));
    core::id::AtomID aid2(rt2->all_bb_atoms().at(4), pose.total_residue());
    std::cout << pose.residue(pose.total_residue()).atom_name(pose.residue(pose.total_residue()).type().all_bb_atoms().at(4)) << std::endl;
    peptide_bb_dof_PHI_term.push_back(core::id::DOF_ID(aid2, core::id::PHI));
    return peptide_bb_dof_PHI_term;
}
std::vector<opt_element> get_peptide_bb_dof_PHI_term_chi2_with_info(core::pose::Pose& pose)
{
    std::vector<opt_element> peptide_bb_dof_PHI_term;
    core::chemical::ResidueTypeOP rt2(new core::chemical::ResidueType(pose.residue(pose.total_residue()).type()));
    core::id::AtomID aid2(rt2->all_bb_atoms().at(4), pose.total_residue());
    std::cout << pose.residue(pose.total_residue()).atom_name(pose.residue(pose.total_residue()).type().all_bb_atoms().at(4)) << std::endl;
    peptide_bb_dof_PHI_term.push_back(

        opt_element(core::id::DOF_ID(aid2, core::id::PHI), std::make_pair(0, 0), pose.residue(pose.total_residue()).type().aa(), "bb phi chi2", 0, pose.total_residue(), 0, 0)

    );
    return peptide_bb_dof_PHI_term;
}
std::vector<core::id::DOF_ID> get_peptide_mc_dof_THETA(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector < core::id::DOF_ID > peptide_mc_dof_THETA;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {

        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 1, atom_end = rt->mainchain_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->mainchain_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(seqpos == 1 && atom == 1)
            {
                continue;
            }
            //std::cout << pose.residue(seqpos).atom_name(rt->mainchain_atoms().at(atom)) << std::endl;
            peptide_mc_dof_THETA.push_back(core::id::DOF_ID(aid, core::id::THETA));
        }
    }
    return peptide_mc_dof_THETA;
}
std::vector<opt_element> get_peptide_mc_dof_THETA_with_info(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector<opt_element> peptide_mc_dof_THETA;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {

        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 1, atom_end = rt->mainchain_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->mainchain_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(seqpos == 1 && atom == 1)
            {
                continue;
            }
            //std::cout << pose.residue(seqpos).atom_name(rt->mainchain_atoms().at(atom)) << std::endl;
            peptide_mc_dof_THETA.push_back(

                opt_element(core::id::DOF_ID(aid, core::id::THETA), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "mc theta", 0, seqpos, 0, 0)

            );
        }
    }
    return peptide_mc_dof_THETA;
}
std::vector<core::id::DOF_ID> get_peptide_mc_dof_D(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector < core::id::DOF_ID > peptide_mc_dof_D;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 1, atom_end = rt->mainchain_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->mainchain_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(seqpos == 1 && atom == 1)
            {
                continue;
            }
            peptide_mc_dof_D.push_back(core::id::DOF_ID(aid, core::id::D));
        }
    }
    return peptide_mc_dof_D;
}
std::vector<opt_element> get_peptide_mc_dof_D_with_info(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector<opt_element> peptide_mc_dof_D;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        core::chemical::ResidueTypeOP rt(new core::chemical::ResidueType(pose.residue(seqpos).type()));
        for(core::Size atom = 1, atom_end = rt->mainchain_atoms().vector().size(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(rt->mainchain_atoms().at(atom), seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            if(seqpos == 1 && atom == 1)
            {
                continue;
            }
            peptide_mc_dof_D.push_back(opt_element(core::id::DOF_ID(aid, core::id::D), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "mc d", 0, seqpos, 0, 0));
        }
    }
    return peptide_mc_dof_D;
}
std::vector<core::id::DOF_ID> get_omega(core::pose::Pose& pose)
{
    std::vector < core::id::DOF_ID > peptide_dof_BB_omega;
    for(core::Size seqpos = 1, seqpos_end = pose.total_residue(); seqpos <= seqpos_end; ++seqpos)
    {
        core::id::TorsionID bb_omg_tid(seqpos, core::id::BB, core::id::omega_torsion);
        peptide_dof_BB_omega.push_back(pose.conformation().dof_id_from_torsion_id(bb_omg_tid));
    }
    peptide_dof_BB_omega.erase(peptide_dof_BB_omega.end());
    return peptide_dof_BB_omega;
}

std::vector<opt_element> get_omega_with_info(core::pose::Pose& pose)
{
    std::vector<opt_element> peptide_dof_BB_omega;
    for(core::Size seqpos = 1, seqpos_end = pose.total_residue(); seqpos <= seqpos_end; ++seqpos)
    {
        core::id::TorsionID bb_omg_tid(seqpos, core::id::BB, core::id::omega_torsion);
        peptide_dof_BB_omega.push_back(
            opt_element(pose.conformation().dof_id_from_torsion_id(bb_omg_tid), std::make_pair(0, 0), pose.residue(seqpos).type().aa(), "bbomega", 0, seqpos, 0, 0));
    }
    peptide_dof_BB_omega.erase(peptide_dof_BB_omega.end());
    return peptide_dof_BB_omega;
}

std::vector<opt_element> get_omega_with_info(core::pose::Pose& pose, std::vector<core::Size> &ind, size_t chain_id, bool remove_last)
{
    std::vector<opt_element> peptide_dof_BB_omega;
    for(core::Size seqpos = 0, seqpos_end = ind.size(); seqpos != seqpos_end; ++seqpos)
    {
        core::id::TorsionID bb_omg_tid(ind[seqpos], core::id::BB, core::id::omega_torsion);
        peptide_dof_BB_omega.push_back(
            opt_element(pose.conformation().dof_id_from_torsion_id(bb_omg_tid), std::make_pair(0, 0), pose.residue(ind[seqpos]).type().aa(), "bbomega", chain_id, ind[seqpos],
                        0, 0));
    }
    if(remove_last)
        peptide_dof_BB_omega.erase(peptide_dof_BB_omega.end());
    return peptide_dof_BB_omega;
}

std::vector<core::id::DOF_ID> get_all_dof_THETA(core::pose::Pose& pose)
{
    std::vector < core::id::DOF_ID > peptide_mc_dof_THETA;
    for(core::Size seqpos = 1, seqpos_end = pose.total_residue(); seqpos <= seqpos_end; ++seqpos)
    {

//        core::chemical::ResidueTypeOP rt( new core::chemical::ResidueType (pose.residue(seqpos).type()) );
//        for(core::Size atom = 1, atom_end = rt->mainchain_atoms().vector().size();
//                atom <= atom_end; ++atom)
        for(core::Size atom = 1, atom_end = pose.residue(seqpos).natoms(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(atom, seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            peptide_mc_dof_THETA.push_back(core::id::DOF_ID(aid, core::id::THETA));
        }
    }
    return peptide_mc_dof_THETA;
}
std::vector<core::id::DOF_ID> get_all_dof_D(core::pose::Pose& pose)
{
    std::vector < core::id::DOF_ID > peptide_mc_dof_D;
    for(core::Size seqpos = 1, seqpos_end = pose.total_residue(); seqpos <= seqpos_end; ++seqpos)
    {
//        core::chemical::ResidueTypeOP rt( new core::chemical::ResidueType (pose.residue(seqpos).type()) );
//        for(core::Size atom = 1, atom_end = rt->mainchain_atoms().vector().size();
//                atom <= atom_end; ++atom)
//        {
//            core::id::AtomID aid(rt->mainchain_atoms().at(atom), seqpos);
        for(core::Size atom = 1, atom_end = pose.residue(seqpos).natoms(); atom <= atom_end; ++atom)
        {
            core::id::AtomID aid(atom, seqpos);
            if(pose.conformation().atom_tree().atom(aid).is_jump())    // skip non-bonded atoms
            {
                continue;
            }
            peptide_mc_dof_D.push_back(core::id::DOF_ID(aid, core::id::D));
        }
    }
    return peptide_mc_dof_D;
}

std::vector<core::id::DOF_ID> get_chi(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector < core::id::DOF_ID > peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            std::cout << "proline! seqpos: " << seqpos << std::endl;
            continue;
        }

        if(pose.residue(seqpos).nchi())
        {
            for(core::Size i = 1, n = pose.residue(seqpos).nchi(); i <= n; ++i)
                peptide_dof_CHI.push_back(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, i)));
        }
    }
    return peptide_dof_CHI;
}

std::vector<core::id::DOF_ID> get_chi1(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector < core::id::DOF_ID > peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            std::cout << "proline! seqpos: " << seqpos << std::endl;
            continue;
        }
        if(pose.residue(seqpos).nchi() >= 1)
        {
            peptide_dof_CHI.push_back(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, 1)));
        }
    }
    return peptide_dof_CHI;
}

std::vector<core::id::DOF_ID> get_chi2(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector < core::id::DOF_ID > peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            std::cout << "proline! seqpos: " << seqpos << std::endl;
            continue;
        }

        if(pose.residue(seqpos).nchi() >= 2)
        {
            peptide_dof_CHI.push_back(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, 2)));
        }
    }
    return peptide_dof_CHI;
}

std::vector<core::id::DOF_ID> get_chi3(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector < core::id::DOF_ID > peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            std::cout << "proline! seqpos: " << seqpos << std::endl;
            continue;
        }

        if(pose.residue(seqpos).nchi() >= 3)
        {
            peptide_dof_CHI.push_back(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, 3)));
        }
    }
    return peptide_dof_CHI;
}

std::vector<core::id::DOF_ID> get_chi4(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector < core::id::DOF_ID > peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            std::cout << "proline! seqpos: " << seqpos << std::endl;
            continue;
        }

        if(pose.residue(seqpos).nchi() == 4)
        {
            peptide_dof_CHI.push_back(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, 4)));
        }
    }
    return peptide_dof_CHI;
}

std::vector<core::id::DOF_ID> get_proline_chi1(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector < core::id::DOF_ID > peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            peptide_dof_CHI.push_back(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, 1)));
        }
    }
    return peptide_dof_CHI;
}
std::vector<opt_element> get_proline_chi1_with_info(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector<opt_element> peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            peptide_dof_CHI.push_back(

                opt_element(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, 1)), std::make_pair(0, 0), pose.residue(seqpos).type().aa(),
                            "proline chi1", 0, seqpos, pose.residue(seqpos).nchi(), 1));
        }
    }
    return peptide_dof_CHI;
}
std::vector<core::id::DOF_ID> get_proline_chi2(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector < core::id::DOF_ID > peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            peptide_dof_CHI.push_back(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, 2)));
        }
    }
    return peptide_dof_CHI;
}
std::vector<opt_element> get_proline_chi2_with_info(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector<opt_element> peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            peptide_dof_CHI.push_back(
                opt_element(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, 2)), std::make_pair(0, 0), pose.residue(seqpos).type().aa(),
                            "proline chi2", 0, seqpos, pose.residue(seqpos).nchi(), 2));
        }
    }
    return peptide_dof_CHI;
}
std::vector<core::id::DOF_ID> get_proline_chi3(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector < core::id::DOF_ID > peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            peptide_dof_CHI.push_back(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, 3)));
        }
    }
    return peptide_dof_CHI;
}

std::vector<opt_element> get_proline_chi3_with_info(core::pose::Pose& pose, core::Size start, core::Size stop)
{
    std::vector<opt_element> peptide_dof_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            peptide_dof_CHI.push_back(

                opt_element(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, 3)), std::make_pair(0, 0), pose.residue(seqpos).type().aa(),
                            "proline chi3", 0, seqpos, pose.residue(seqpos).nchi(), 3));
        }
    }
    return peptide_dof_CHI;
}

std::vector<std::tuple<core::id::DOF_ID, core::chemical::AA, core::Size, core::Size, core::Size> > get_peptide_all_chi_dof_with_aatype(core::pose::Pose& pose,
        core::Size start,
        core::Size stop)
{
    std::vector<std::tuple<core::id::DOF_ID, core::chemical::AA, core::Size, core::Size, core::Size> > peptide_dof_and_aatype_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if(pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr)
        {
            std::cout << "proline! seqpos: " << seqpos << std::endl;
            continue;
        }
        if(pose.residue(seqpos).nchi() >= 1)
        {
            for(core::Size j = 1; j <= pose.residue(seqpos).nchi(); j++)
            {
                peptide_dof_and_aatype_CHI.push_back(
                    std::make_tuple(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, j)), pose.residue(seqpos).type().aa(),
                                    pose.residue(seqpos).nchi(), j, seqpos));
            }
        }
    }
    return peptide_dof_and_aatype_CHI;
}

std::vector<opt_element> get_peptide_all_chi_dof_with_info(core::pose::Pose& pose, core::Size start, core::Size stop, bool include_proline)
{
    std::vector<opt_element> peptide_dof_and_aatype_CHI;
    for(core::Size seqpos = start, seqpos_end = stop; seqpos <= seqpos_end; ++seqpos)
    {
        if((pose.residue(seqpos).type().aa() == core::chemical::aa_pro || pose.residue(seqpos).type().aa() == core::chemical::aa_dpr) && !include_proline)
        {
            std::cout << "proline! seqpos: " << seqpos << std::endl;
            continue;
        }
        if(pose.residue(seqpos).nchi() >= 1)
        {
            for(core::Size j = 1; j <= pose.residue(seqpos).nchi(); j++)
            {
                peptide_dof_and_aatype_CHI.push_back(opt_element(
                        pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(seqpos, core::id::CHI, j)),
                        std::make_pair(-numeric::NumericTraits<core::Real>::pi(), numeric::NumericTraits<core::Real>::pi()),
                        pose.residue(seqpos).type().aa(),
                        "chi", 0,
                        seqpos, pose.residue(seqpos).nchi(), j));
            }
        }
    }
    return peptide_dof_and_aatype_CHI;
}

std::vector<opt_element> get_peptide_all_chi_dof_with_info(core::pose::Pose& pose, std::vector<core::Size> &ind, size_t chain_id, bool include_proline)
{
    std::vector<opt_element> peptide_dof_and_aatype_CHI;
    for(core::Size seqpos = 0, seqpos_end = ind.size(); seqpos != seqpos_end; ++seqpos)
    {
        if((pose.residue(ind[seqpos]).type().aa() == core::chemical::aa_pro || pose.residue(ind[seqpos]).type().aa() == core::chemical::aa_dpr) && !include_proline)
        {
            std::cout << "proline! seqpos: " << ind[seqpos] << std::endl;
            continue;
        }
        if(pose.residue(ind[seqpos]).nchi() >= 1)
        {
            for(core::Size j = 1; j <= pose.residue(ind[seqpos]).nchi(); j++)
            {
                peptide_dof_and_aatype_CHI.push_back(
                    opt_element(pose.conformation().dof_id_from_torsion_id(core::id::TorsionID(ind[seqpos], core::id::CHI, j)),
                                std::make_pair(-numeric::NumericTraits<core::Real>::pi(), numeric::NumericTraits<core::Real>::pi()), pose.residue(ind[seqpos]).type().aa(), "chi",
                                chain_id, ind[seqpos], pose.residue(ind[seqpos]).nchi(), j));
            }
        }
    }
    return peptide_dof_and_aatype_CHI;
}

void insert_to_opt_vector(std::vector<opt_element> &opt_vect_with_info,
                          core::pose::Pose &pose,
                          std::string arguments,
                          pepsgo::ranges &range)
{
    if(arguments.find("phipsi") != std::string::npos)
    {
        std::vector<opt_element> peptide_mc_p_info = get_phi_psi_with_info(pose);

        size_t i1 = 0, i2 = 0, zero = 0;
        i1 = opt_vect_with_info.size(), i2 = opt_vect_with_info.size() + peptide_mc_p_info.size() - 1;
        range.phipsi =
            (peptide_mc_p_info.size()) ? std::make_tuple(true, i1, i2) : std::make_tuple(false, zero, zero);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_mc_p_info.begin(), peptide_mc_p_info.end());
        // peptide_mc_p.pop_back(); // NH2 to C-term
    }
    if(arguments.find("omega") != std::string::npos)
    {
        std::vector<opt_element> peptide_mc_o_info = get_omega_with_info(pose);

        size_t i1 = 0, i2 = 0, zero = 0;
        i1 = opt_vect_with_info.size(), i2 = opt_vect_with_info.size() + peptide_mc_o_info.size() - 1;
        range.omega =
            (peptide_mc_o_info.size()) ? std::make_tuple(true, i1, i2) : std::make_tuple(false, zero, zero);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_mc_o_info.begin(), peptide_mc_o_info.end());
    }
    if(arguments.find("allchi") != std::string::npos)
    {
        std::vector<opt_element> all_chi_info =
            get_peptide_all_chi_dof_with_info(pose, 1, pose.total_residue(), true);

        size_t i1 = 0, i2 = 0, zero = 0;
        i1 = opt_vect_with_info.size(), i2 = opt_vect_with_info.size() + all_chi_info.size() - 1;
        range.chi = (all_chi_info.size()) ? std::make_tuple(true, i1, i2) : std::make_tuple(false, zero, zero);
        opt_vect_with_info.insert(opt_vect_with_info.end(), all_chi_info.begin(), all_chi_info.end());
    }
    if(arguments.find("allchiexceptpro") != std::string::npos)
    {
        std::vector<opt_element> all_chi_info =
            get_peptide_all_chi_dof_with_info(pose, 1, pose.total_residue(), false);

        size_t i1 = 0, i2 = 0, zero = 0;
        i1 = opt_vect_with_info.size(), i2 = opt_vect_with_info.size() + all_chi_info.size() - 1;
        range.chi = (all_chi_info.size()) ? std::make_tuple(true, i1, i2) : std::make_tuple(false, zero, zero);
        opt_vect_with_info.insert(opt_vect_with_info.end(), all_chi_info.begin(), all_chi_info.end());
    }
    if(arguments.find("proline") != std::string::npos)
    {
        std::vector<opt_element> peptide_proline_chi1 = get_proline_chi1_with_info(pose, 1, pose.total_residue());
        std::vector<opt_element> peptide_proline_chi2 = get_proline_chi2_with_info(pose, 1, pose.total_residue());
        std::vector<opt_element> peptide_proline_chi3 = get_proline_chi3_with_info(pose, 1, pose.total_residue());

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_proline_chi1.begin(), peptide_proline_chi1.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_proline_chi2.begin(), peptide_proline_chi2.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_proline_chi3.begin(), peptide_proline_chi3.end());
    }
    if(arguments.find("sctheta12") != std::string::npos)
    {
        std::vector<opt_element> peptide_sc_t_h =
            get_peptide_sc_dof_THETA_with_info(pose, 1, pose.total_residue(), 1);
        std::vector<opt_element> peptide_sc_t_a =
            get_peptide_sc_dof_THETA_with_info(pose, 1, pose.total_residue(), 2);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_sc_t_h.begin(), peptide_sc_t_h.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_sc_t_a.begin(), peptide_sc_t_a.end());
    }
    if(arguments.find("bbtheta12") != std::string::npos)
    {
        std::vector<opt_element> peptide_bb_t_h =
            get_peptide_bb_dof_THETA_without_mc_with_info(pose, 1, pose.total_residue(), 1);
        std::vector<opt_element> peptide_bb_t_a =
            get_peptide_bb_dof_THETA_without_mc_with_info(pose, 1, pose.total_residue(), 2);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_bb_t_h.begin(), peptide_bb_t_h.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_bb_t_a.begin(), peptide_bb_t_a.end());
    }
    if(arguments.find("mctheta") != std::string::npos)
    {
        std::vector<opt_element> peptide_mc_t = get_peptide_mc_dof_THETA_with_info(pose, 1, pose.total_residue());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_mc_t.begin(), peptide_mc_t.end());

        // peptide_mc_t.erase(peptide_mc_t.begin()); // ACETYLATED_NTERMINUS_VARIANT
    }

    if(arguments.find("scd12") != std::string::npos)
    {
        std::vector<opt_element> peptide_sc_d_h = get_peptide_sc_dof_D_with_info(pose, 1, pose.total_residue(), 1);
        std::vector<opt_element> peptide_sc_d_a = get_peptide_sc_dof_D_with_info(pose, 1, pose.total_residue(), 2);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_sc_d_h.begin(), peptide_sc_d_h.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_sc_d_a.begin(), peptide_sc_d_a.end());
    }
    if(arguments.find("bbd12") != std::string::npos)
    {
        std::vector<opt_element> peptide_bb_d_h =
            get_peptide_bb_dof_D_without_mc_with_info(pose, 1, pose.total_residue(), 1);
        std::vector<opt_element> peptide_bb_d_a =
            get_peptide_bb_dof_D_without_mc_with_info(pose, 1, pose.total_residue(), 2);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_bb_d_h.begin(), peptide_bb_d_h.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_bb_d_a.begin(), peptide_bb_d_a.end());
    }
    if(arguments.find("mcd") != std::string::npos)
    {
        std::vector<opt_element> peptide_mc_d = get_peptide_mc_dof_D_with_info(pose, 1, pose.total_residue());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_mc_d.begin(), peptide_mc_d.end());
    }

    if(arguments.find("scp12") != std::string::npos)
    {
        std::vector<opt_element> peptide_sc_p_h =
            get_peptide_sc_dof_PHI_without_CHI_with_info(pose, 1, pose.total_residue(), 1);
        std::vector<opt_element> peptide_sc_p_a =
            get_peptide_sc_dof_PHI_without_CHI_with_info(pose, 1, pose.total_residue(), 2);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_sc_p_h.begin(), peptide_sc_p_h.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_sc_p_a.begin(), peptide_sc_p_a.end());
    }
    if(arguments.find("bbp12") != std::string::npos)
    {
        std::vector<opt_element> peptide_bb_p_h =
            get_peptide_bb_dof_PHI_without_mc_without_term_with_info(pose, 1, pose.total_residue(), 1);
        std::vector<opt_element> peptide_bb_p_a =
            get_peptide_bb_dof_PHI_without_mc_without_term_with_info(pose, 1, pose.total_residue(), 2);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_bb_p_h.begin(), peptide_bb_p_h.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_bb_p_a.begin(), peptide_bb_p_a.end());
    }

    if(arguments.find("frp12") != std::string::npos)
    {
        std::vector<opt_element> peptide_fr_p_h = get_peptide_bb_dof_PHI_first_term_with_info(pose, 1);
        std::vector<opt_element> peptide_fr_p_a = get_peptide_bb_dof_PHI_first_term_with_info(pose, 2);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_fr_p_h.begin(), peptide_fr_p_h.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_fr_p_a.begin(), peptide_fr_p_a.end());
    }
    if(arguments.find("lrp12") != std::string::npos)
    {
        std::vector<opt_element> peptide_lr_p_h = get_peptide_bb_dof_PHI_last_term_with_info(pose, 1);
        std::vector<opt_element> peptide_lr_p_a = get_peptide_bb_dof_PHI_last_term_with_info(pose, 2);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_lr_p_h.begin(), peptide_lr_p_h.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_lr_p_a.begin(), peptide_lr_p_a.end());
    }

    if(arguments.find("12chi") != std::string::npos)
    {
        std::vector<opt_element> peptide_2_c = get_peptide_bb_dof_PHI_term_chi2_with_info(pose);
        std::vector<opt_element> peptide_1_c = get_peptide_bb_dof_PHI_term_chi1_with_info(pose);

        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_2_c.begin(), peptide_2_c.end());
        opt_vect_with_info.insert(opt_vect_with_info.end(), peptide_1_c.begin(), peptide_1_c.end());
    }
}

}
