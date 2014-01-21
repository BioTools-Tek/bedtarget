#include "targeter.h"


void Targeter::targetPromoters(short margin_ds){
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        uint txStart = gh->txStart;

        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        //Region before
        uint txStart_before = txStart - margin_ds;

        cout << chrom << '\t'
             << txStart_before << '\t' << txStart << '\t'
             << gene_name << "|PROMOTER_" << margin_ds << endl;
    }
}


//Score is defined in grabexons

void Targeter::printExtras(QChar &score1, QChar &score2, QChar myway, short frame){
    if(this->scores){
        //Append Quality
        QList<QChar> &keys = checks;
        cout << '\t' << flush;

        if(keys.contains(score1))  cout << score1.toLatin1() << flush;
        else cout << 'w' << flush;

        if(keys.contains(score2))  cout << score2.toUpper().toLatin1() << flush;
        else cout << 'W' << flush;
    }

    if (this->direction){
        cout << '\t' << myway.toLatin1() << flush;
    }

    // -99 is null
    if (this->frames){
        cout << '\t' << frame << flush;
    }
    cout << endl;
}



void Targeter::targetExons()
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        for(int j=0; j < gh->exons.length(); j++)
        {
            ExonHolder *ex = gh->exons.at(j);
            //Dont print UTR
            if (ex->utr==-1)
            {
                cout << chrom << '\t'
                     << ex->start << '\t' << ex->stop << '\t'
                     << gene_name << '|'
                     << Exn << gh->determineExonNumber(ex->exon) << flush; printExtras(gh->score_start,gh->score_endl, (gh->direction?'+':'-'), ex->frame);
            }
        }
    }
}


void Targeter::targetUTRExons()
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        for(int j=0; j < gh->exons.length(); j ++)
        {
            ExonHolder *ex = gh->exons.at(j);

            // Print UTR ONLY
            // No need to worry about splice sites

            if(ex->utr!=-1){

                cout << chrom << '\t'
                     << ex->start << '\t' << ex->stop << '\t'
                     << gene_name << '|'
                     << Exn << gh->determineExonNumber(ex->exon) << flush;

                if(ex->utr==5) cout << '_' << UTR5;
                else cout << '_' << UTR3;

                cout << flush; printExtras(gh->score_start,gh->score_endl, (gh->direction?'+':'-'), ex->frame);
            }
        }
    }
}


void Targeter::targetUTRRegions()
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        cout << chrom << '\t'
             << gh->txStart << '\t' << gh->cdsStart
             << '\t' << gene_name  << '|'
             << (gh->direction?UTR5:UTR3) << flush;printExtras(gh->score_start,gh->score_endl);

        cout << chrom << '\t'
             << gh->cdsStop << '\t' << gh->txStop
             << '\t' << gene_name  << '|'
             << (gh->direction?UTR5:UTR3) << flush; printExtras(gh->score_start,gh->score_endl);

    }
}



void Targeter::targetIntergenic(bool codingOnly)
{
    QList<GeneHolder*> &genes = this->genelist;
    uint last_genepos = 0; // Intergenic position holder
    QString last_gene = "";

    int gen_size = genes.size();
    for(int i=0; i < gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);

        uint tmp_start = codingOnly?gh->cdsStart:gh->txStart;
        uint tmp_stop = codingOnly?gh->cdsStop:gh->txStop;

        if((i==0) || (last_gene.split("-ISO")[0] == gh->gene_name.split("-ISO")[0])){
            last_genepos = tmp_stop;
            last_gene = gh->gene_name;
            continue;
        }

        cout << gh->chrom.toUtf8().data() << '\t' << flush;

        //Isoforms can be fickle
        if(last_genepos > tmp_start){
            cout << tmp_start << '\t' << last_genepos
                 << '\t' << gh->gene_name.toUtf8().data() << '-' << last_gene.toUtf8().data() << flush;
        } else {
            cout << last_genepos << '\t' << tmp_start
                 << '\t' << last_gene.toUtf8().data() << '-' << gh->gene_name.toUtf8().data() << flush;
        }
        cout << '|' << Ing << flush; printExtras(gh->score_start,gh->score_endl);

        last_genepos = tmp_stop;
        last_gene = gh->gene_name;
    }
}



void Targeter::targetGenes(bool codingOnly)
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for(int i=0; i < gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);

        uint tmp_start = codingOnly?gh->cdsStart:gh->txStart;
        uint tmp_stop = codingOnly?gh->cdsStop:gh->txStop;

        cout << gh->chrom.toUtf8().data() << '\t'
             << tmp_start << '\t' << tmp_stop  << '\t'
             << gh->gene_name.toUtf8().data()
             << flush; printExtras(gh->score_start,gh->score_endl);
    }
}


void Targeter::targetSpliceOnly(Splice *ss, bool no_utr)
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        for(int j=0; j < gh->exons.length(); j ++)
        {
            ExonHolder *ex = gh->exons.at(j);

            // If dont print utr, and exon is not utr, skip.. wat
            if (no_utr && ex->utr!=-1) continue;
 //         if (ex->utr==-1)
 //         {

            //Print SPLICE ONLY
            //First exon, only splice on cdsStop side
            //Last exon, only splice on cdsStart side

            uint splice1_start = ex->start;
            uint splice1_end = ex->start;

            uint splice2_start = ex->stop;
            uint splice2_end = ex->stop;

            // Bedtarget recieves the exons in +'ve order
            // so the last exon is the same in any given direction

            // txStart cdsStart  first-exon
            //+       ______5'UTR OOOOOOOOOOO 3'E__A____5'S  --->
            //- <---  ______3'UTR OOOOOOOOOOO 5'E__D____3'S

            //               last exon     cdsStop    txStop
            //+         5'S__D___5'E OOOOOOOOOOO 3'UTR_____ --->
            //-  <---   3'E__A___3'S OOOOOOOOOOO 5'UTR_____

            //the placing of splice_sites is actually independent of direction.
            // it depends only on first and last index

            if (!(ex->start==gh->cdsStart)) splice1_start -= ss->splice_site;
            if (!(ex->stop==gh->cdsStop)) splice2_end += ss->splice_site;

            //Every now and then there are UTR_exons so we should ignore any
            //             txStart             cdsStart
            // +      5'UTR_______DonSP ooooooo|OOOOOOOOO AccSP_____
            //splice sites at cdsStart or cdsStop, but not any inbetween ones

            // The order of printing splice sites should still be splice1 followed
            // by splice2. Just the naming of donor and acceptor change.

            //if direction is +ve: donor_sites = splice1, accept_sites = splice2
            //if direction is -ve: accept_sites= splice1,  donor_sites = splice2

            QString utr_mention="";
            if(ex->utr==5) utr_mention = QString('_').append(UTR5);
            else if (ex->utr==3) utr_mention = QString('_').append(UTR3);


            if(gh->direction){
                if(ss->donor_sites){
                    cout << chrom << '\t'
                         << splice1_start << '\t' << splice1_end << '\t'
                         << gene_name << '|'
                         << Exn << gh->determineExonNumber(ex->exon) << utr_mention.toUtf8().data() << '_' << DonSpl << flush; printExtras(gh->score_start,gh->score_endl, (gh->direction?'+':'-'), ex->frame);
                }
                if(ss->acceptor_sites){
                    cout << chrom << '\t'
                         << splice2_start << '\t' << splice2_end << '\t'
                         << gene_name << '|'
                         << Exn << gh->determineExonNumber(ex->exon) << utr_mention.toUtf8().data() << '_' << AccSpl << flush; printExtras(gh->score_start,gh->score_endl, (gh->direction?'+':'-'), ex->frame);
                }
            }
            else{
                if(ss->acceptor_sites){
                    cout << chrom << '\t'
                         << splice1_start << '\t' << splice1_end << '\t'
                         << gene_name << '|'
                         << Exn << gh->determineExonNumber(ex->exon) << utr_mention.toUtf8().data() << '_' << AccSpl << flush; printExtras(gh->score_start,gh->score_endl, (gh->direction?'+':'-'), ex->frame);
                }
                if(ss->donor_sites){
                    cout << chrom << '\t'
                         << splice2_start << '\t' << splice2_end << '\t'
                         << gene_name << '|'
                         << Exn << gh->determineExonNumber(ex->exon) << utr_mention.toUtf8().data() << '_' << DonSpl << flush; printExtras(gh->score_start,gh->score_endl, (gh->direction?'+':'-'), ex->frame);
                }
            }
        }
    }
}


void Targeter::targetIntrons(Splice *ss)
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {

        uint last_exonpos = 0;

        GeneHolder *gh = genes.at(i);
        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        for(int j=0; j < gh->exons.length(); j ++)
        {
            ExonHolder *ex = gh->exons.at(j);

            uint splice1_start = ex->start;
            uint splice2_end = ex->stop;

            if (!(ex->start==gh->cdsStart)) splice1_start -= ss->splice_site;
            if (!(ex->stop==gh->cdsStop)) splice2_end += ss->splice_site;

            if(j==0) {
                last_exonpos = splice2_end;
                continue;
            }

            cout << chrom << '\t'
                 << last_exonpos << '\t' << splice1_start << '\t'
                 << gene_name << '|'
                 << Itr << (gh->determineExonNumber(ex->exon)-1) << flush; printExtras(gh->score_start,gh->score_endl, (gh->direction?'+':'-'), ex->frame);
        }
    }
}






