#ifndef ARGS_H
#define ARGS_H

#define VERSION "2.4_20170521"

#include <QString>
#include <iostream>
#include <stdio.h>
#include "splice.h"
#include <getopt.h>

#define assert(x){ cerr << x << endl; cerr << endl; exit(-1);}

using namespace std;

// REQ ARGS
#define ARG_GENES "genes"
#define ARG_EXONS "exons"
#define ARG_AUTOSOMES "autosomes"

//Opts
#define ARG_HELP "help"

#define ARG_DATABASE "database"
#define ARG_TABLE "table"
#define ARG_LOCAL "local"

#define ARG_NONCODING "include-noncoding"
#define ARG_INTRONS "introns"
#define ARG_INTERGENIC "intergenic"
#define ARG_SCORES "scores"
#define ARG_FRAMES "frames"
#define ARG_DIRECTION "direction"
#define ARG_UTR "utr"
#define ARG_UTRONLY "utr-only"
#define ARG_PROMOTERS "promoters"
#define ARG_SPLICE "splice"
#define ARG_SPLICEONLY "splice-only"
#define ARG_SPLICENOACCEPT "splice-no-acceptors"
#define ARG_SPLICENODONOR "splice-no-donors"
#define ARG_REFSEQGENE "refseq-gene-names"

#define ARG_STANDARD "standard"

// Standard is the equivalent of
//exons,scores,frames,directions,utr,splice 5


static struct option long_options[] =
{
{ARG_DATABASE,        required_argument, 0,   'b'},
{ARG_HELP,            no_argument,       0,   'h'},
{ARG_TABLE,           required_argument, 0,   't'},
{ARG_LOCAL,           no_argument,       0,   'l'},
{ARG_UTR,             no_argument,       0,   'u'},
{ARG_UTRONLY,         no_argument,       0,   'U'},
{ARG_NONCODING,       no_argument,       0,   'N'},
{ARG_INTRONS,         no_argument,       0,   'n'},
{ARG_INTERGENIC,      no_argument,       0,   'r'},
{ARG_PROMOTERS,       required_argument, 0,   'p'},
{ARG_SPLICE,          required_argument, 0,   's'},
{ARG_SPLICEONLY,      required_argument, 0,   'S'},
{ARG_SPLICENOACCEPT,  no_argument,       0,   'A'},
{ARG_SPLICENODONOR,   no_argument,       0,   'D'},
{ARG_GENES,           no_argument,       0,   'g'},
{ARG_EXONS,           no_argument,       0,   'e'},
{ARG_SCORES,          no_argument,       0,   'k'},
{ARG_FRAMES,          no_argument,       0,   'f'},
{ARG_DIRECTION,       no_argument,       0,   'd'},
{ARG_REFSEQGENE,      no_argument,       0,   'q'},
{ARG_STANDARD,        no_argument,       0,   'x'},
{0,                   0,                 0,    0 }
};

class Options{
public:
    QString database, table;
    bool local, codingonly, introns, intergenic, scores;
    bool exonframes, direction, unique_isoform_names;
    short promoter_margin_downstream, promoter_margin_upstream, utrlevel;
    Splice *ss;

    Options() :
        database(""), table(""),
        local(false), codingonly(false),
        introns(false), intergenic(false), scores(false), exonframes(false),
        direction(false), unique_isoform_names(false),
        promoter_margin_downstream(-1), promoter_margin_upstream(-1),
        utrlevel(-1), ss(0)
    {}

};





struct Args {
    QStringList region_list;
    bool genes, exons;

    Options opt_arg;



    QList<int> region_num_check(QString nums)
    {
        QStringList r1_r2 = nums.split('-');


        int first=-1,second=-1;
        bool r1b;

        first = r1_r2.first().split("hr").last().toInt(&r1b);
        if (!r1b) assert("Please give a valid number for the start of the range, e.g. chr1-22, or chr1:12345-");

        if (r1_r2.size()==2){
            second = r1_r2.last().toInt(&r1b);
            if(!r1b) assert("Please give a valid number for the end of the range, e.g. chr1-22, or chr1:2-123456");
            if (second <= first) assert("End of range cannot be smaller than start of range!");
        }
        else if (r1_r2.size()>2) assert("Please give a single range, no extra hyphens '-'");

        if (second==-1) second=22;

        return QList<int> ({first, second});
    }

    void usage(char * prog)
    {
        fprintf(stderr,
                "version %s\n"
                "Usage:\t%s <region> (--genes or --exons) [OPTIONS]\n\n"
                "  or\t%s --help\n\t    for more information\n\n\n",
                VERSION, prog, prog);
        exit(-1);

    }


    void help(char *prog){
        fprintf(stderr,
                "version %s\n"
                "Usage:\t%s <region> (--genes or --exons) [OPTIONS]"
                "\n\n"
                "REGION FORMAT:\n"
                "    chrN[:reg1-[reg2]]\n"
                " or chrN-[M]    (where N<M and M <= 22)\n"
                " or %s   (equiv. \"chr1-22\")\n"
                " (multiple ranges can be specified seperated by \"region1,region2,...\")\n\n"
                "MAP OPTIONS:\n"
                "  --%s\tor -%s [hg18|hg19|etc]\n"
                "  --%s\tor -%s [refGene]\n"
                "  --%s\tor -%s\n"
                "\n"
                "GENE OPTIONS:\n"
                "  --%s\tor -%s,\n\t(Default switch: shortcut for\n"
                "         --splice 5 --utr --scores --frames\n"
                "         --direction --refseq-gene-names)\n\n"
                "  --%s\tor -%s\n"
                "  --%s\tor -%s, (REMOVES ALL OTHER FLAGS)\n"
                "  --%s\tor -%s\n"
                "  --%s\tor -%s <downstream_margin> <upstream_margin>\n"
                "  --%s\tor -%s\n"
                "  --%s\tor -%s\n"
                "  --%s\tor -%s <margin> \n"
                "  --%s\tor -%s, (REMOVES ALL OTHER NON-SPLICE FLAGS)\n"
                "       --%s\tor -%s\n"
                "       --%s\tor -%s\n"
                "  --%s\tor -%s\n"
                "  --%s\tor -%s\n"
                "  --%s\tor -%s\n"
                "  --%s\tor -%s\n",
                VERSION, prog,
                ARG_AUTOSOMES,
                ARG_DATABASE, "b",
                ARG_TABLE, "t", ARG_LOCAL, "l",

                ARG_STANDARD, "x",
                ARG_UTR, "u", ARG_UTRONLY, "U",
                ARG_NONCODING, "N",
                ARG_PROMOTERS, "p",
                ARG_INTRONS, "n",
                ARG_INTERGENIC, "r",
                ARG_SPLICE, "s", ARG_SPLICEONLY, "S", ARG_SPLICENOACCEPT, "A", ARG_SPLICENODONOR, "D",
                ARG_SCORES, "k",
                ARG_FRAMES, "f",
                ARG_DIRECTION, "d",
                ARG_REFSEQGENE, "q"
                );
        exit(-1);
    }



    void handleRegionArg(char * reg, char * prog)
    {
        QString region = reg;

        // Catch help text here
        if (region.startsWith("--" ARG_HELP)) help(prog);

        QStringList original_region_list = region.split(',');

        for (QString region : original_region_list)
        {

            if (region==ARG_AUTOSOMES) region="chr1-22";

            if (region.startsWith("chr"))
            {
                if (region.contains(':')) // Single region with bp range
                {
                    QString nums = region.split("hr").last().split(':').last();
                    if (nums.contains('-')){
                        region_num_check(nums);
                    }
                    else assert("Please give a valid region: chr1:12345-, or chr2:145-19901, etc");
                }
                else { //No ':', single chrom (entire) or chrom range

                    region = region.split("hr").last();

                    if (region.contains('-')){
                        QList<int> small_big = region_num_check(region);

                        int small_chrom = small_big[0];
                        int big_chrom = small_big[1];

                        for (int i=small_chrom; i <= big_chrom; i++)
                            region_list.append(QString("chr").append(QString::number(i)));

                        //                        delete small_big;
                        continue;
                    }
                    else { // Just a single chromosome
                        region_list.append(QString("chr").append(region));
                        continue;
                    }
                }
            }
            else assert("Invalid region, start with 'chr' at least");

            // Made it this far, append
            region_list.append(region);

        }
    }



    Args(int argc, char ** argv) : genes(false), exons(false)
    {
        bool utrtoo = false;
        bool utronly = false;
        int splice_amount=-1;
        bool splice_only = false, splice_acc = true, splice_don = true;
        bool splice_at_all=false;

        if (argc<2) usage(argv[0]);

        // == Required Arguments ==
        handleRegionArg(argv[1], argv[0]);

        if (argc<3) usage(argv[0]);


        // --genes or --exons
        QString second = argv[2];
        if (second.startsWith("--")){
            if (second == "--" ARG_EXONS) exons = true;
            else if (second == "--" ARG_GENES) genes = true;
            else if (second == "--" ARG_HELP) help(argv[0]);
            else usage(argv[0]);
        }
        else usage(argv[0]);


        extern char *optarg;
        int ch;

        int option_index = 3;

        while((ch=getopt_long(argc, argv, "hxb:t:luUNnrp:s:SADkfdq", long_options, &option_index))!=-1){

            switch(ch){
            case 'h':
                help(argv[0]); exit(-1); break;
            case 'b':
                opt_arg.database = optarg;  break;
            case 't':
                opt_arg.table = optarg;     break;
            case 'l':
                opt_arg.local = true;       break;
            case 'u':
                utrtoo = true;      break;
            case 'U':
                utronly=true;       break;
            case 'N':
                opt_arg.codingonly=false;   break;
            case 'n':
                opt_arg.introns = true;     break;
            case 'r':
                opt_arg.intergenic = true;  break;
            case 's':
                splice_amount = QString(optarg).toInt(&splice_at_all);
                if (!splice_at_all){
                    cerr << "Give a splice margin value. Essential sites are 2, and highly conserved are 5 bp." << endl;
                    exit(-1);
                }
                break;

            case 'p':
                bool conv1, conv2;
                opt_arg.promoter_margin_downstream = QString(argv[optind-1]).toShort(&conv1);
                opt_arg.promoter_margin_upstream = QString(argv[optind]).toShort(&conv2);

                if (not (conv1 || conv2)){
                    cerr << "Give downsream and upstream promoter values. Default is 5000 and 0." << endl;
                    exit(-1);
                }
                break;

            case 'S':
                splice_only = true;  break;
            case 'A':
                splice_acc = false; break;
            case 'D':
                splice_don = false; break;
            case 'g':
                genes = true;       break;
            case 'e':
                exons = true;       break;
            case 'k':
                opt_arg.scores = true;      break;
            case 'f':
                opt_arg.exonframes = true;      break;
            case 'd':
                opt_arg.direction = true;   break;
            case 'q':
                opt_arg.unique_isoform_names = true;
                break;

            case 'x':
                //                scores,frames,directions,utr,splice 5
//                exons = true;
                opt_arg.scores= true;
                opt_arg.exonframes = true;
                opt_arg.direction = true;

                opt_arg.unique_isoform_names = true;
                utrtoo = true;
                splice_at_all = true;
                splice_amount = 5; //Essential
                break;

            case ':':
                cerr << "missing opt" << endl;
                exit(-1);
                break;

            default:
                usage(argv[0]);
                break;
            }
        }

        //Check tables
        if (opt_arg.database==""){
            cerr << "Assuming database: hg19" << endl;
            opt_arg.database = "hg19";
        }

        if (opt_arg.table==""){
            cerr << "Assuming table: refGene" << endl;
            opt_arg.table = "refGene";
        }

        //Check opts
        if (genes && exons) assert("Either exons or genes, one or the other");

        if(!genes && !exons){ // Nothing specified
            if(!opt_arg.introns && !opt_arg.intergenic && opt_arg.promoter_margin_downstream==-1){
                usage(argv[0]);
                //              assert("nothing specified");
            }
        }

        //Splice
        if(splice_at_all && !exons) assert("--exons must be specified if --splice is used");
        if ((splice_acc || splice_don) && splice_amount==0) assert("--splice must be specified too");
        //UTR
        if (genes){
            if (utronly) opt_arg.utrlevel = 2;
            else if (utrtoo) opt_arg.utrlevel = 1;
            else opt_arg.utrlevel = 0;
        }
        else if (exons){
            if (utronly) opt_arg.utrlevel = 2;
            else if (utrtoo) opt_arg.utrlevel = 1;
            else opt_arg.utrlevel = 0;
        }
        opt_arg.ss = new Splice(splice_amount,splice_don, splice_acc, splice_only);

        //        delete optarg;

//        debug();
    }


    void debug(){
        cerr << "regions="; for (QString region : region_list) cerr << region.toUtf8().data() << ",";
        cerr << endl;
        cerr << "genes=" << genes << ", exons=" << exons << endl;
        cerr << "\nOpts:" << endl;
        cerr << "db=" << opt_arg.database.toUtf8().data()
             << ", table=" << opt_arg.table.toUtf8().data() << ", local=" << opt_arg.local << endl;
        cerr << "codingonly=" << opt_arg.codingonly << ", introns=" << opt_arg.introns
             << ", intergenic=" << opt_arg.intergenic << ", scores=" << opt_arg.scores
             << ", exonframes=" << opt_arg.exonframes << ", direct=" << opt_arg.direction
             << ", unique_names=" << opt_arg.unique_isoform_names << endl;
        cerr << "promoterU=" << opt_arg.promoter_margin_upstream << ", D=" << opt_arg.promoter_margin_downstream << endl;
        cerr << "UTR=" << opt_arg.utrlevel << endl;
        cerr << "Splice[acc don only marg]=" << opt_arg.ss->acceptor_sites << " "
             << opt_arg.ss->donor_sites << " " << opt_arg.ss->splice_only
             << " " << opt_arg.ss->splice_site << endl;
//        exit(-1);
    }



};





#endif // ARGS_H



