#include "grabexons.h"


QList<QString> GrabExons::retrieveFromMYSQL()
{
    //Make command string
    QString command;
    if (local) command = QString("mysql -uroot -peggwax -A --execute=\"SELECT * FROM ");
    else command = QString("mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A --execute=\"SELECT * FROM ");

    command.append(database).append('.').append(table)
            .append(" WHERE ")
            .append("chrom='").append(chrom)
            .append('\'');

    if(reg1!=-1 && reg2!=-1) command.append(" AND cdsStart BETWEEN ").append(QString::number(reg1)).append(" AND ").append(QString::number(reg2)).append(' ');
    else command.append(' ');
    command.append("ORDER BY cdsStart;\"");

    //Process and get output
    QProcess *qp = new QProcess; //Not a child process.
    qp->start(command);
    qp->waitForFinished(90000);

    QString all_data = qp->readAllStandardOutput();
    QList<QString> tmp = all_data.split('\n');

    return tmp;
}


GeneHolder *GrabExons::parseExons(const QString &lines)
{
    QStringList data = lines.split('\t');
    // chrom=2, strand=3, exonCount=8, starts=9,ends=10, name2=12, exonFrames=15, startstat=13, endstat=14
    //THESE FIELDS ARE NEEDED AS DESCRIBED HERE https://lists.soe.ucsc.edu/pipermail/genome/2006-November/012218.html
    QStringList startsR = data.at(9).split(','), endsR = data.at(10).split(','), framesR = data.at(15).split(',');

    GeneHolder *gh = new GeneHolder(data.at(12).trimmed()); //Name

    if (gh->gene_name=="name2"){
        gh->chrom="0";
        return gh;
    }

    //Coding region
    gh->cdsStart = data.at(6).toInt();
    gh->cdsStop = data.at(7).toInt();

    gh->txStart = data.at(4).toInt();
    gh->txStop  = data.at(5).toInt();

    gh->exonCount = data.at(8).toInt();
    gh->direction = ((data.at(3).trimmed()).at(0)=='+');
    gh->chrom = data.at(2).trimmed();

    //Parse reading frames
    for (int p=0; p< framesR.size(); p++){
        gh->reading_frames[p] = QString(framesR[p]).toShort();
    }


    if(startsR.length() != endsR.length()){
        cerr << "starts and ends length dont match" << endl;
        exit(-1);
    }



    QList<ExonHolder*> exons;

    for(int i=0; i < startsR.length(); i++)
    {
        uint start = startsR.at(i).toInt();
        uint end = endsR.at(i).toInt();

        if ((start==0) && (end==0)){
            continue;// Random... why is this even here?
        }

        short utr = -1;
        short frame = gh->reading_frames[i];

        //Scenarios to check:
        //1. Completely before codingStart --> 5' UTR
        //2. Completely after codingEnd -->  3'UTR
        //3. Start before codingStart, end after codingStart --> split
        //4. Start before codingEnd, end after codingEnd --> split

        //1.
        if (end <= gh->cdsStart){
            utr = (gh->direction)?5:3;
            exons.append(new ExonHolder(start,end, (i+1), utr, frame));
        }
        //2.
        else if(start >= gh->cdsStop){
            utr= (gh->direction)?3:5;
            exons.append(new ExonHolder(start,end, (i+1), utr, frame));
        }
        else if(utr==-1){          //UTR not found, search partial
            //       5'UTR | Exn
            //3. +---oooooo|OOOOOOOO---+
            if ( (start < gh->cdsStart) && (end > gh->cdsStart))
            {
                exons.append(new ExonHolder(start,gh->cdsStart, (i+1), (gh->direction)?5:3, frame));  //BeforeUTR
                exons.append(new ExonHolder(gh->cdsStart,end,(i+1), -1, frame));                      //After Exon

            }
            //4. +----OOOOOOOO|ooo-----+
            else if ( (start < gh->cdsStop) && (end > gh->cdsStop))
            {
                exons.append(new ExonHolder(start,gh->cdsStop, (i+1), -1, frame));                 //BeforeExon
                exons.append(new ExonHolder(gh->cdsStop,end,(i+1), ((gh->direction)?3:5), frame)); //AfterUTR;
            }
            //Partial not found --> Coding Exon!
            else {
                exons.append(new ExonHolder(start,end, (i+1), -1, frame));
            }
        }
        else {
            cerr << "NOTHING SHOULD BE HERE" << endl;
        }
    }

    //Scores
    gh->score_start = data.at(13).trimmed().at(0);
    gh->score_endl = data.at(14).trimmed().at(0);
//    cerr << gh->score_start.toLatin1() << ',' << gh->score_endl.toLatin1() << endl;

    gh->exons = exons;                //Not a reference

    //Determine unique exon numbers
    QList<uint> unique_exons;
    for (int ue=0; ue < exons.length(); ue ++){
        uint num = exons[ue]->exon;
        if (!unique_exons.contains(num)) unique_exons.append(num);
    }

    //Assumes exons are in order
    gh->exon_numbers = unique_exons;


    //Check for isoforms
    int isonum = isoformcheck(gh);
    if (isonum > 0) gh->gene_name.append("-ISOF").append(QString::number(isonum));

    return gh;
}


uint GrabExons::isoformcheck(GeneHolder *gh)
{
    QString &chrom = gh->chrom;
    QString &name = gh->gene_name;

    if (!chromemap.contains(chrom)) {
        QMap<QString, uint> tmp;
        chromemap[chrom] = tmp;
    }

    if(!chromemap[chrom].contains(name)){
        return (chromemap[chrom][name] = 0); // Count how many times this gene has been inserted
    }
    else return (++chromemap[chrom][name]);
}
