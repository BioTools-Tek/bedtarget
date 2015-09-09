#ifndef TARGETER_H
#define TARGETER_H

#include "grabexons.h"
#include "Args.h"
#include <QList>


class Targeter
{
private:
    typedef struct {
        uint pos1;
        uint pos2;
    } pair;
public:
    QList<GeneHolder*> genelist;
    bool scores, frames, direction;

    Targeter(QList<GeneHolder*> &genelist, Args &args)
    {

        this->genelist = genelist;
        this->scores = args.opt_arg.scores;
        this->frames = args.opt_arg.exonframes;
        this->direction = args.opt_arg.direction;

        if(args.opt_arg.promoter_margin_downstream!=-1) targetPromoters(\
                    args.opt_arg.promoter_margin_upstream,\
                    args.opt_arg.promoter_margin_downstream); // targets genes only

        if(args.opt_arg.intergenic) targetIntergenic(args.opt_arg.codingonly);

        if(args.opt_arg.introns) targetIntrons(args.opt_arg.ss);


        if(args.genes){
            if(args.opt_arg.utrlevel == 0) targetGenes(args.opt_arg.codingonly);
            else if(args.opt_arg.utrlevel == 1) {
                targetGenes(args.opt_arg.codingonly);
                targetUTRRegions();
            }
            else if(args.opt_arg.utrlevel == 2) targetUTRRegions();
        }

        if(args.exons){

            // Include splice sites if non-zero splice margins set
            if(args.opt_arg.ss->splice_site!=0) {
                targetSpliceOnly(args.opt_arg.ss, !(args.opt_arg.utrlevel==1));
                // last argument
                //no_utr = false
                //no_utr = true
            }


            // Also...
            if(!args.opt_arg.ss->splice_only){
                if(args.opt_arg.utrlevel == 0) targetExons();
                else if(args.opt_arg.utrlevel == 1){
                    targetExons();
                    targetUTRExons();
                }
            }
            else if(args.opt_arg.utrlevel == 2) targetUTRExons();

        }
    }

    //Target Jobs
    //Key:
    //    ~ -- coding in progress
    //    ' -- coding done, testing needed
    //    " -- testing
    //    ¬ -- bugs, needs fixing
    //    ^ -- working, complete

    //amount to subtract from middle exons
    //amount to subtract from middle exons
    void targetIntergenic(bool codingOnly);   //'
    void targetGenes(bool codingOnly);        //'
    void targetExons();                       //'
    void targetIntrons(Splice *ss);           //'
    void targetUTRExons();                    //'
    void targetUTRRegions();                  //'
    void targetSpliceOnly(Splice *ss, bool no_utr=true);        //'
    void targetSpliceUTR(Splice *ss);         //'
    void targetPromoters(short margin_us, short margin_ds);

    //Other jobs
    void sortedList();
    void printList();
    void printExtras(GeneHolder *&gh, short frame=-99);

};

#endif // TARGETER_H
