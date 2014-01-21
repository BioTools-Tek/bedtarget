#include "targeter.h"
#include <string>

int main(int argc, char *argv[])
{
    Args argh(argc, argv);

    GrabExons ex(argh.database, argh.region,
                 argh.table, argh.local);

    QList<GeneHolder*> &genes = ex.genes;
    cerr << argh.region.toStdString() << "\t#Exons: " << ex.genes.length() << endl;

    if (genes.size()==0){
        cerr << "[Error] No genes found - Check your connection" << endl;
        exit(-1);
    }

    Targeter(genes, argh);

}
