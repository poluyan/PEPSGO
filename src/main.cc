#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

int main(int argc, char *argv[]){
    std::cout << "Start..." << std::endl;
    devel::init(argc, argv);
    std::cout << "Hello wold!" << std::endl;
    core::pose::Pose peptide;
    std::string peptide_seq = "AAAAA";
    core::pose::make_pose_from_sequence(
	peptide, peptide_seq, 
	*core::chemical::ChemicalManager::get_instance()->residue_type_set( 
	core::chemical::FA_STANDARD));

    //core::scoring::ScoreFunctionOP scorefn = core::scoring::ScoreFunctionFactory::create_score_function(core::scoring::REF_2015);
    core::scoring::ScoreFunctionOP scorefn = core::scoring::get_score_function();

    (*scorefn)(peptide);
    std::cout << peptide.energies().total_energy() << std::endl;
    peptide.dump_pdb("ala5.pdb");
}