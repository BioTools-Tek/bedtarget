#ifndef ARGS_H
#define ARGS_H

#define VERSION "1.9c"

#include <QString>
#include <iostream>
#include <stdio.h>
#include "splice.h"
#include <getopt.h>

#define assert(x){ cerr << x << endl; cerr << endl; exit(-1);}

using namespace std;

// ARGS
#define ARG_DATABASE "database"
#define ARG_TABLE "table"
#define ARG_LOCAL "local"

#define ARG_GENES "genes"
#define ARG_EXONS "exons"
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

// Standard is the equivalent of\
exons,scores,frames,directions,utr,splice 5


static struct option long_options[] =
{
    {ARG_DATABASE,        required_argument, 0,   'b'},
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
    bool local, codingonly, genes, exons, introns, intergenic, scores, exonframes, direction, unique_isoform_names;
    short promoter_margin_downstream, promoter_margin_upstream, utrlevel;
    Splice *ss;
};


struct Args {
    QString region;
    Options opt_arg;


    void usage(char *prog){
        fprintf(stderr,
        "version %s\n"
        "Usage: %s <chrN[:reg1-reg2]> [MAP] [GENE]\n"
        "\n"
        "MAP OPTIONS:\n"
        "  --%s [hg18|hg19|etc]\n"
        "  --%s [refGene]\n"
        "  --%s\n"
        "\n"
        "GENE OPTIONS:\n"
        "  --%s\n"
        "  --%s\n"
        "  --%s\n"
        "  --%s <downstream_margin> <upstream_margin>\n"
        "  --%s\n"
        "  --%s\n"
        "  --%s <margin> \n"
        "  --%s\n"
        "   --%s\n"
        "   --%s\n"
        "  --%s || --%s\n"
        "  --%s\n"
        "  --%s\n"
        "  --%s\n"
        "  --%s\n"
        "  --%s\n",
        VERSION,
        prog,
        ARG_DATABASE,  ARG_TABLE, ARG_LOCAL,
        ARG_UTR, ARG_UTRONLY,
        ARG_NONCODING,
        ARG_PROMOTERS,
        ARG_INTRONS,
        ARG_INTERGENIC,
        ARG_SPLICE, ARG_SPLICEONLY, ARG_SPLICENOACCEPT, ARG_SPLICENODONOR,
        ARG_GENES, ARG_EXONS,
        ARG_SCORES,
        ARG_FRAMES,
        ARG_DIRECTION,
        ARG_STANDARD,
        ARG_REFSEQGENE);
        exit(-1);
    }



    Args(int argc, char ** argv){

        bool utrtoo = false;
        bool utronly = false;
        int splice_amount=0;
        bool splice_only = false, splice_acc = true, splice_don = true;


        if (argc<2) usage(argv[0]);
        region = argv[1];

        extern char *optarg;
        int ch;

        int option_index = 0;

        while((ch=getopt_long(argc, argv, "x:b:t:luUNnp:rs:S:ADgekfdq", long_options, &option_index))!=-1){

            switch(ch){
            case'b':
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
                splice_amount = QString(optarg).toInt();
                if (splice_amount==0){
                    cerr << "Give a splice value. Default is 5." << endl;
                    exit(-1);
                }
                break;

            case 'p':
                bool conv1, conv2;
                opt_arg.promoter_margin_downstream = QString(argv[option_index+1]).toShort(&conv1);
                opt_arg.promoter_margin_upstream = QString(argv[option_index+2]).toShort(&conv2);

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
                opt_arg.genes = true;       break;
            case 'e':
                opt_arg.exons = true;       break;
            case 'k':
                opt_arg.scores = true;      break;
            case 'f':
                opt_arg.exonframes = true;      break;
            case 'd':
                opt_arg.direction = true;   break;
            case 'q':
                opt_arg.unique_isoform_names = true;
                break;
            default:
                usage(argv[0]);
                break;
            }
        }

        //Check opts
        if (opt_arg.genes && opt_arg.exons) assert("Either exons or genes, one or the other");

        if(!opt_arg.genes && !opt_arg.exons){ // Nothing specified
            if(!opt_arg.introns && !opt_arg.intergenic && opt_arg.promoter_margin_downstream==-1){
                usage(argv[0]);
    //              assert("nothing specified");
            }
        }

        //Splice
        if(splice_amount!=0 && !opt_arg.exons) assert("--exons must be specified if --splice is used");
        if ((splice_acc || splice_don) && splice_amount==0) assert("--splice must be specified too");
        //UTR
        if (opt_arg.genes){
            if (utronly) opt_arg.utrlevel = 2;
            else if (utrtoo) opt_arg.utrlevel = 1;
            else opt_arg.utrlevel = 0;
        }
        if (opt_arg.exons){
            if (utronly) opt_arg.utrlevel = 2;
            else if (utrtoo) opt_arg.utrlevel = 1;
            else opt_arg.utrlevel = 0;
        }
        opt_arg.ss = new Splice(splice_amount,splice_don, splice_acc, splice_only);
    }
};





#endif // ARGS_H



