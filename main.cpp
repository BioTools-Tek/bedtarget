#include "targeter.h"
#include <string>


bool GeneHolder::print_direction;
bool GeneHolder::print_scores;
bool GeneHolder::print_uniquenames;
bool GeneHolder::print_frames;



int main(int argc, char *argv[])
{
    Args argh(argc, argv);

    GeneHolder::print_direction = argh.opt_arg.direction;
    GeneHolder::print_scores = argh.opt_arg.scores;
    GeneHolder::print_uniquenames = argh.opt_arg.unique_isoform_names;
    GeneHolder::print_frames = argh.opt_arg.exonframes;

    GeneHolder::printHeader();

    //Handle Data
    for (QString region : argh.region_list)
    {
        GrabExons ex(argh.opt_arg.database, region,
                     argh.opt_arg.table, argh.opt_arg.local);

        QList<GeneHolder*> &genes = ex.genes;
        cerr << region.toStdString() << "\t#Exons: " << ex.genes.length() << endl;

        if (genes.size()==0){
            cerr << "[Error] No genes found - Check your connection" << endl;
            exit(-1);
        }

        Targeter(genes, argh);
    }
}
