#-------------------------------------------------
#
# Project created by QtCreator 2013-06-12T13:35:38
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = bedtarget2
CONFIG   += console
CONFIG   += c++11
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    grabexons.cpp \
    targeter.cpp

HEADERS += \
    grabexons.h \
    targeter.h \
    Args.h \
    splice.h \
    mysql.h
