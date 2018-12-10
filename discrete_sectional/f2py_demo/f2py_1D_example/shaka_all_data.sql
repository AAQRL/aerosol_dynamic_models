DROP TABLE IF EXISTS all_data CASCADE;
CREATE TABLE all_data (
id          SERIAL UNIQUE PRIMARY KEY,
date        timestamp with time zone,
username    varchar(12),
hostname    varchar(40),
ipnbr       varchar(15),
os          varchar(20),     -- output of 'uname' to get platform
sysinfo     varchar(2500),   -- text blob for output of sysinfo
runtime     float,           -- wallclock sec to complete
signal      varchar(80),     -- signal file (path/filename)
theta       float,           -- zero, order phase correction
nmr_methd   varchar(10),     -- FDM, RRT or DFT
Nsig        int,             -- signal length to be used
wmino       int,             -- window size min value
wmaxo       int,             -- window size max value
par         varchar(80),     -- linelist output file
threshhold  float,           -- threshhold for par above
ReSp        varchar(80),     -- spectra output files
ImSp        varchar(80),     -- spectra output files
AbsSp       varchar(80),     -- spectra output files
rho         float,           -- basis density
Nb0         int,             -- basis size
Nbc         int,             -- coarse basis
Nsp         int,             -- plotting pnts
Gamm        float,           -- smoothing (5d-2)
cheat       float,           -- used in FDM, see discriptions.
cheatmore   char(1),         -- used in FDM, see discriptions.
ros         float            -- RRT regularization parameter
);
