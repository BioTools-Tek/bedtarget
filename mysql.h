#ifndef MYSQL_H
#define MYSQL_H

#include <QString>
#include <QProcess>
#include <iostream>

using namespace std;

class MySQL {
public:
    static QString retrieveMYSQLResult(QString query){
        QString command = "mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A --execute=\"";
        command.append(query).append("\"");

        //Process and get output
        QProcess *qp = new QProcess(); //Not a child process.

    #ifdef DEBUG
        cerr << command.toUtf8().data() << endl;
    #endif
        qp->start(command);
        qp->waitForFinished(900000000);
        QString all_data = qp->readAllStandardOutput();

    #ifdef DEBUG
       cerr << all_data.toUtf8().data() << endl;
    #endif

       delete qp;
       return all_data;
    }
};

#endif // MYSQL_H
