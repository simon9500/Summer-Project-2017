{\rtf1\ansi\ansicpg1252\cocoartf1504\cocoasubrtf830
{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red255\green255\blue255;}
{\*\expandedcolortbl;;\csgray\c0;\csgray\c100000;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 INCLUDED FILES:\
\
findSN.py (finds exact location of SN, no uncertainties)\
sample_dt.py (samples start times from the distribution and produces a histogram)\
sampleplotSN.py (samples start times, applies them to dt values and plots the circles)\
SNfunctions.py (a set of functions used in my python programs)\
astropy (I added the astropy library in the same directory as a hack to get it imported correctly)\
t1.out\
t1thresh.out\
t2.out\
t2thresh.out\
t3.out\
t3thresh.out\
t4.out\
t4thresh.out\
\
The t1.out (etc.) files are the digitised graphs of the time distributions.  The t1thresh.out (etc.) files are the digitised energy threshold/integrated rate graphs.\
\
REQUIRED LIBRARIES:\
(these are the versions I have that work on my computer)\
\
numpy 
\f1\fs22 \cf2 \cb3 \CocoaLigature0 1.13.1\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf2 astropy 
\f1\fs22 2.0.1
\f0\fs24 \
matplotlib 
\f1\fs22 2.0.2\

\f0\fs24 ROOT 
\f1\fs22 6.08/06 
\f0\fs24 (I installed it using {\field{\*\fldinst{HYPERLINK "https://alexpearce.me/2016/02/root-on-os-x-el-capitan/"}}{\fldrslt https://alexpearce.me/2016/02/root-on-os-x-el-capitan/}})\
\
HOW TO RUN:\
\
While in the directory where the files are type the following in the terminal:\
\

\f1\fs22 python findSN.py\
python sample_dt.py t1.out t1thresh.out\
python sampleplotSN.py t1.out t1thresh.out\
\

\f0\fs24 To test findSN.py with different SN locations or dt values you will need to edit the source code to change the variables.  An Aitoff graph will be saved at the end showing the circles produced and the estimated SN location will also be displayed in terminal output.\
\
Similarly, to run sampleplotSN.py with different event numbers and for different numbers of start times edit the source code.  An Aitoff histogram will be produced at the end.
\f1\fs22 \

\f0\fs24 \
\

\f1\fs22 \
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 \cb1 \CocoaLigature1 \
\
}