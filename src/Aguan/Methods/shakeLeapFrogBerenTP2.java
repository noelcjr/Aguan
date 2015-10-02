/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Aguan.Methods;
import Aguan.TheMatrix.TheMatrix;
import Aguan.files.log;
import Aguan.files.restart;
import java.io.*;

/**
 *
 * @author noel
 */
class gbl {                           /* globals structure */
    double dt;                          /* time step */
    double dtsq;                        /* timeStepSquared */
    double sigma;                       /* for Leonard-Jones */
    double epslon;                      /* for Leonard-Jones */
    double rcut;                        /* cut off for PE */
    double rcutsq;                      /* " squared */
    double rcutsqbig;                   /* rcut+2*bondLength squared */
    double tolerance;                   /* for constraint loop */
    double sigsq;                       /* leonard-jones */
    double boxinv;                      /* 1/box */
    double box;                         /* length of box */
    double wnc, wc;                     /* unconstrained, constrained virial */
    double temperature;                 /* for temperature scaling */
    double tt;                          /* coupling coef for temp .1ps */
    double pt;                          /* coupling coef for press .1ps */
    double desiredpressure;             /* for pressure scaling */
    double isothermalbeta;              /* for pressure scaling */
    double j2kcal;                      /* convert joules to kcal/mole */
    double pscalefactor;                /* mu (scale rx,ry,rz, volume) */
    double tscalefactor;                /* lambda (scale vx,vy,vz) */
    THERMO w;                           /* total virial */
    THERMO boxvolume;                   /* box cubed */
    THERMO ke;                          /* kinetic */
    THERMO pe;                          /* potential energy */
    THERMO energy;                      /* pe+ke */
    THERMO calcpressure;                /* calculated pressure */
    THERMO compressibilityfactor;       /* pv/nkt */
    THERMO calctemperature;             /* temperature: calc */
    THERMO hbonds;                      /* hydrogen bond stats */
    int[] hbondbins;                    /* histogram of hbonds/molecule */
    int paircnt;                        /* for paircnt averaging */
    double bindelta;                    /* width of radial pair bin */
    int vacfptr;                        /* for veloc auto correl */
    int evalcnt;                        /* how often to do pair dist, etc */
    int ifaceflag;                      /* set if in user iface routine */
    int savecnt;                        /* for file saving, round robin */
    int cnt, maxcnt;                    /* iteration count, and max count */
    int MAXHBONDBINS;
    public gbl(int maxhbondbins){
        MAXHBONDBINS = maxhbondbins;
        hbondbins = new int[MAXHBONDBINS];
    }
}
class THERMO{
    /* for thermodynamic variables */
    double cur;             /* current value */
    double min;             /* min value found */
    double max;             /* max value found */
    double tot;             /* sum values */
    double rms;             /* sum values squared */
    double cnt;             /* total values */
}

public class shakeLeapFrogBerenTP2 {
 public int HY1 = 0;                   /* indexes into structure array */
 public int OXY = 1;
 public int HY2 = 2;
 public String[] atomname = {"HY1","OXY","HY2"};
 public int[] nxtitab = {HY1,HY1,HY1,OXY,OXY,HY2,HY2,HY2,OXY};
 public int[] nxtjtab = {HY1,OXY,HY2,HY1,HY2,HY1,OXY,HY2,OXY};
 public int[] binseltab = {0,2,1,2,1,1,2,1,2};
 //#define UNIX 0              /* set if unix (ie gcc) */
 public int DODUMP = 1;            /* ascii to stdout */
 public int DODUMPTSTEP = 0;       /* display (ascii) KE,PE,etc every time step */
 public int DODUMPTSTEPLOW = 0;    /* ascii, rx,ry,rz,vx,vy,vz per timestep */
 public int DODUMPVERBOSE = 0;     /* ascii debug */
 public int DODUMPERRORS = 0;      /* ascii debug */
 public int FORCE_NVE = 0;         /* 1=overide Pressure/Temp scaling */
 public double DESIRED_TEMPERATURE = 300.0;    /* kelvin  (only for NPT) */
 public double DESIRED_PRESSURE = 1.00e5;      /* pascal (1e5 = 1bar) (only for NPT) */
 public int SPC_E = 0;             /* 1=use spc/e value for charge */
 public int NPT2NVE = 0;           /* 1= Continue Run, but clear arrays, thermo vals */
    // Works with local variables only
 //   public double[] rxinit;
 public int chkflags;
 public double piconst,boltzmanconstant;    /* pi, etc... */
 public double avogadroconst, elementarycharge, permittivityconst, deg2rad;
 public char contname;
 public int INBUFSIZE = 300;
 public char[] holdbuf = new char[INBUFSIZE];
 public char[] Konsbuf1 = new char[INBUFSIZE];
 public char[] Konsbuf2 = new char[INBUFSIZE];
 public int TEMPERATURECONSTRAINT = 3;
 public int NUMMOLECULES = (6 * 6 * 6);        /* total # of water molecules */
 public int NUMCONSTRAINTS = 3;                /* rigid SPC */
 public int NUMATOMS = 3;                      /* h-o-h */
 public int NUMCOORD = 3;                      /* x,y,z */
 public int NUMFREEDOM = (NUMMOLECULES * (NUMCOORD * NUMATOMS - NUMCONSTRAINTS) - TEMPERATURECONSTRAINT);
 public int NUMPARTICLES = (NUMMOLECULES * NUMATOMS);
 public double GRAMS2KGRAMS = 0.001;
 public double KILO = 1000.0;
 public double ANGSTROM = 1.0e-10;             /* 1angstrom = 1e-10 meters */
 public double PICOSEC = 1.0e12;
 public double PASCAL2ATM = (1.0/1.013e5);
 public double PASCAL2BAR = (1.0/1.000e5);
 public double JOULE2CAL = 0.2389;

 public int MAXKCNT = 100;                 /* constraint loop abort, see ltolerence */
 public int MAXBIN = 500;                  /* num bins for pair distribution */
 public int LOC_WHITE = 15;                /* graphics colors */
 public int LOC_RED = 14;
 public int LOC_GREEN = 13;
 public int LOC_BLUE = 12;
 public int LOC_YELLOW = 11;
 public int TARGETMOLECULE = 122;              /* watch this mol., for hbonds, etc */
 public double POTENTIALHYDROBOND = -2.80;     /* was 2.25, -2.0 */
 public int EVALVACF = 5;                      /* evaluate VACF every 5 iterations */
 public int WAITAVG = 5;
 public int EVALCNTMAX = 10;                   /* eval pair dist every 10 iter */
 public int MAXEPDIST = 550;                   /* max bins */
 public int NUMPAIRS = 3;                      /* o-o o-h h-h */
 public int NUMVACFBINS = 400;
 public int ACFNUMREF = 4;
 public int VACFSPACING = NUMVACFBINS / ACFNUMREF;  /* 400/4 = 100 */
 public int MAXHBONDBINS = 10;
 public double EPDISTMIN = 7.0;
 public double EPDISTMULT = 50.0;
 public double MASSHYDROGEN = 1.00797;
 public double MASSOXYGEN = 15.9994;
 public double BIGNUMBER = 1.0e99;
 //struct _pair {              /* radial pair distribution structure */
 public int[][] pairs_bin = new int[MAXBIN][NUMPAIRS];         /* 0:oo, 1:oh 2:hh */
 public double[][] pairs_gn = new double[MAXBIN][NUMPAIRS];    /* after smoothing filter */
 //};struct _pair pairs[MAXBIN];
 public int[] energypairdist = new int[MAXEPDIST];      /* energy pair distribution array */
 //struct _acfvel{             /* for Velocity auto-correl function */
 public double acfvel_vx[][] = new double[ACFNUMREF][NUMMOLECULES];   /* hold reference vx0 */
 public double acfvel_vy[][] = new double[ACFNUMREF][NUMMOLECULES];   /* hold reference vy0 */
 public double acfvel_vz[][] = new double[ACFNUMREF][NUMMOLECULES];   /* hold reference vz0 */
 public int[] acfvel_validflag = new int[ACFNUMREF];
//};struct _acfvel acfvel[ACFNUMREF];

//struct _vacf {              /* for velocity auto-corr. function */
 public double[] vacf_total = new double[NUMVACFBINS];           /* sum of correl values */
 public int[] vacf_cnt = new int[NUMVACFBINS];                /* total count */
//};struct _vacf vacf[NUMVACFBINS];

 private gbl globalStructure;
 private THERMO tmr;
 //Water
 //atom[] atoms = new atom[3];                    /* 0:hydrogen, 1:oxygen, 2:hydrogen */
 public double[][] water_dsq = new double[NUMCONSTRAINTS][3];     /* the squared bond length constraints */
 public double[][] water_rx = new double[NUMMOLECULES][3];
 public double[][] water_ry = new double[NUMMOLECULES][3];
 public double[][] water_rz = new double[NUMMOLECULES][3];     /* current positions  r(t) */
 public double[][] water_vx = new double[NUMMOLECULES][3];
 public double[][] water_vy = new double[NUMMOLECULES][3];
 public double[][] water_vz = new double[NUMMOLECULES][3];     /* current velocity v(t-.5dt) */
 public double[][] water_fx = new double[NUMMOLECULES][3];
 public double[][] water_fy = new double[NUMMOLECULES][3];
 public double[][] water_fz = new double[NUMMOLECULES][3];     /* forces f(t) */
 public double[][] water_mass = new double[NUMMOLECULES][3];
 public double[][] water_charge = new double[NUMMOLECULES][3];
 public int[] water_kcolor = new int[NUMMOLECULES];
 //double lastx,lasty;  /* for screen update only */
 //END WATER
 private log LOG;
public shakeLeapFrogBerenTP2(TheMatrix TM, restart RR){
        chkflags = 0;
        contname = 'a';
        initUniverseConstants();
        System.out.println("Size of init array = "+rxinit.length);
        System.out.println("contname = "+contname);
        globalStructure = new gbl(10);
        tmr = new THERMO();
        LOG = new log(TM.logfile);
        initmolecules(contname,globalStructure);
        md(globalStructure);
}

private void initmolecules(char contname, gbl globStruc){
        double chargehydrogen, chargeoxygen, masshydrogen, massoxygen;
        double bondlength, bondangle, virtualbondlength;
        int i,j,kk;

        /* basic values */
    globStruc.j2kcal = JOULE2CAL * avogadroconst /(KILO * (double)NUMMOLECULES); /* kcal/mole of water molecules */
    globStruc.dt = 5.0e-16;                           /* seconds (was 2e-15) */
    globStruc.maxcnt = 800;                         /*80000 = 40ps */
    globStruc.rcut = 8.5 * ANGSTROM;                  /* was 8.5 */
    globStruc.tolerance = 1.0e-6;                     /* was e-6 e-5 old:e-4 */
    globStruc.temperature = DESIRED_TEMPERATURE;      /* kelvin */
    globStruc.desiredpressure = DESIRED_PRESSURE;     /* pascal */
    globStruc.tt = 0.100e-12;              /* .1ps [berendsen 1984] */
    globStruc.pt = 0.100e-12;              /* .1ps [berendsen 1984] */
    globStruc.isothermalbeta = 2.0e-10;    /* [berendsen 1984] was 4.9e-10 */
    globStruc.box = 18.62 * ANGSTROM;      /* stillinger   density = 1gm/cm3 */

    /* SPC water (berendsen j.phys.chem 1987 91,6269)
     * and from prevost et al  mol.phys 1990 vol.71, #3 587-603
     *  a=.37122, b=.3428, q=.41 (or .4238)
     */
    globStruc.sigma = 3.16555789 * ANGSTROM;         /* meters (water)  (3.4argon) */
    globStruc.epslon = .65 * KILO / avogadroconst;   /* joules (water) =1.08e-21 (1.65e-21argon) */
    bondlength = 1.0 * ANGSTROM;               /* SPC water  H-O */
    bondangle = 109.47;                        /* SPC water  degrees  H-O-H */
    chargehydrogen = .41 * elementarycharge;   /*  coulomb */
    if(SPC_E == 1){
    chargehydrogen = .4238 * elementarycharge; /* coulomb */
    }
    chargeoxygen = -2.0 * chargehydrogen;    /* coulomb */

    /* other variables */
    globStruc.bindelta = .02 * ANGSTROM;      /* width of radialPair bin */
    globStruc.dtsq = globStruc.dt * globStruc.dt;                 /* secondsSq */
    virtualbondlength = 2.0 * bondlength * Math.sin( deg2rad * bondangle / 2.0);  /* H-H */
    masshydrogen = (MASSHYDROGEN / avogadroconst) * GRAMS2KGRAMS;  /* Kgrams */
    massoxygen = (MASSOXYGEN / avogadroconst) * GRAMS2KGRAMS;  /* Kgrams */
    globStruc.rcutsq = globStruc.rcut * globStruc.rcut;
    globStruc.rcutsqbig = (globStruc.rcut + 2.0*bondlength) * (globStruc.rcut + 2.0*bondlength);
    globStruc.sigsq = globStruc.sigma * globStruc.sigma;
    globStruc.savecnt = 0;
    globStruc.cnt = 0;
    globStruc.wnc = 0.0;
    globStruc.wc = 0.0;
    globStruc.ifaceflag = 0;

    initlow(1, globStruc);

    kk = 0;
    for (i=0; i < NUMMOLECULES; i++){   /* init the water structure */
        water_dsq[i][0] = bondlength * bondlength;
        water_dsq[i][1] = bondlength * bondlength;
        water_dsq[i][2] = virtualbondlength * virtualbondlength;
        //    water[i].kcolor = LOC_BLUE;
        /* printf("i:%d dsq0:%lg dsq1:%lg dsq2:%lg\n",i,water[i].dsq[0],water[i].dsq[1],water[i].dsq[2]); */
        for (j=0; j < NUMATOMS; j++){
            //water[i].atoms[j].lastx = 0.0;
            //water[i].atoms[j].lasty = 0.0;
            water_rx[i][j] = rxinit[kk++];
            water_ry[i][j] = rxinit[kk++];
            water_rz[i][j] = rxinit[kk++];
            water_vx[i][j] = rxinit[kk++];
            water_vy[i][j] = rxinit[kk++];
            water_vz[i][j] = rxinit[kk++];
            water_fx[i][j] = 0;
            water_fy[i][j] = 0;
            water_fz[i][j] = 0;
            if (j == OXY){
                water_mass[i][j] = massoxygen;
                water_charge[i][j] = chargeoxygen;
            }
            else{
                water_mass[i][j] = masshydrogen;
                water_charge[i][j] = chargehydrogen;
            }
            LOG.lout.printf("%d %d %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f \n",i,j,water_rx[i][j],water_ry[i][j],water_rz[i][j],water_vx[i][j],water_vy[i][j],water_vz[i][j],water_mass[i][j],water_charge[i][j]);
        }
    }

    //if(s){    /* continue from previous run (if any) */
    //    restorestate(s);
    // #if NPT2NVE
    //     gbl.boxvolume.cur = gbl.boxvolume.tot/gbl.boxvolume.cnt;    /* fixup */
    //     gbl.box = pow( gbl.boxvolume.cur , (1/3.0));          /* for pressure */
    //     initlow(0);     /* continue from NPT, clear arrays, thermo values */
    //#endif
    //}

    globStruc.boxinv = 1.0/globStruc.box;

    //sprintf(Konsbuf1,"mass(H:%lg O:%lg) charge(H:%lg O:%lg) ",masshydrogen,massoxygen,chargehydrogen,chargeoxygen);  /* kludge, hold for disptherm routine */
    //sprintf(Konsbuf2,"bondlength:%lg bondangle:%lg virtualbondlen:%lg  ",bondlength,bondangle,virtualbondlength);    /* ditto */
 }
private void md(gbl glb){
     double t, kselfdiff, hpe, hdensity;
     kselfdiff = 0;
     for(t=glb.box*glb.dt; glb.cnt < glb.maxcnt; t += glb.dt, glb.cnt++){
         if((glb.cnt % 1) == 0)
             savestate(glb);
         glb.evalcnt++;
         if(glb.evalcnt == EVALVACF)
             createacfvelbins(glb);
         if(glb.evalcnt >= EVALCNTMAX){
             kselfdiff = calcselfdiff(t,glb);
             calcpairdist(glb);
             glb.evalcnt = 0;
             createacfvelbins(glb);
         }
         hpe = force(glb);
         evalscalefactors(hpe,glb);
         mover(glb);               /* leap-frog with constraints (calc ke, wc) */

        glb.box = glb.box * glb.pscalefactor;  /* pressure scale, [berendsen et al, 1984, step 9p] */
        glb.boxinv = 1.0/glb.box;
        upthermo(glb.boxvolume, glb.box * glb.box * glb.box, glb.cnt);
        upthermo(glb.energy, glb.pe.cur + glb.ke.cur, glb.cnt);
        hdensity = (double)NUMMOLECULES * (2.0 * MASSHYDROGEN + MASSOXYGEN) /(glb.boxvolume.cur * avogadroconst * 1.0e6);
    if(DODUMPTSTEP == 1){
        LOG.lout.printf("%d: PE:%8.6g  KE:%8.6g  E:%8.6gKcal T:%8.6g  ",glb.cnt,glb.pe.cur*glb.j2kcal,glb.ke.cur*glb.j2kcal,glb.energy.cur*glb.j2kcal, glb.calctemperature.cur);
        //LOG.lout.rtmess(holdbuf,2,1);
        LOG.lout.printf(" Pressure:%8.6g PV/NKT:%8.6g Virial:%8.6g  ",glb.calcpressure.cur,glb.compressibilityfactor.cur,glb.w.cur);
        //LOG.lout.rtmess(holdbuf,2,2);
        LOG.lout.printf(" SelfDiff:%8.6g  Box:%8.6g  Density:%8.6g   ",kselfdiff, glb.box, hdensity);
        //prtmess(holdbuf,2,3);
        if(DODUMP == 1){
            LOG.lout.printf("%d\n",glb.cnt);
        }
     }
    //    dumper(t);
     }

}
private double[][] mat = new double[3][3];
private double[][] inv = new double[3][3];
private void mover(gbl glb){
    //ATOM *aptr;
    int i,kcnt;

    double rx1, ry1, rz1, rx2, ry2, rz2, rx3, ry3, rz3;  /* current posit */
    double px1, py1, pz1, px2, py2, pz2, px3, py3, pz3;  /* next posit */
    double rx12, ry12, rz12, rx23, ry23, rz23, rx31, ry31, rz31; /* current relative posit */
    double px12, py12, pz12, px23, py23, pz23, px31, py31, pz31; /* current next posit */
    double om1, om2, om3;            /* 1/mass */
    double r12sq, r23sq, r31sq;      /* bond constraints */
    double r12r23, r23r31, r31r12;   /* scalar products */
    double p12sq, p23sq, p31sq;
    double p12r12, p12r23, p12r31;
    double p23r12, p23r23, p23r31;
    double p31r12, p31r23, p31r31;
    double l12, l23, l31, l12new, l23new, l31new;
    double const1, const2, const3;
    double quad1, quad2, quad3, vec1, vec2, vec3, determ;
    double q111, q122, q133, q112, q123, q131;
    double q211, q222, q233, q212, q223, q231;
    double q311, q322, q333, q312, q323, q331;
    double ltolerance = glb.dtsq * glb.tolerance;
    double kenertmp = 1.0/glb.dt;
    double vx,vy,vz;
    double hke;

    hke = 0.0;
    glb.wc = 0.0;
    /* for atoms: 0:hydrogen, 1:oxygen, 2:hydrogen */
    for (i=0; i < NUMMOLECULES; i++){
        r12sq = water_dsq[i][0];       //water[i].dsq[0];      /* bond length constraints */
        r23sq = water_dsq[i][1];
        r31sq = water_dsq[i][2];

        /* leap-frog algorithm */
        //aptr = &water[i].atoms[HY1];   /* 1st hydrogen */
        om1 = 1.0/(water_mass[i][HY1]);
        vx = water_vx[i][HY1] + water_fx[i][HY1] * om1 * glb.dt;       /* v(t+.5dt) */
        vy = water_vy[i][HY1] + water_fy[i][HY1] * om1 * glb.dt;
        vz = water_vz[i][HY1] + water_fz[i][HY1] * om1 * glb.dt;
        vx *= glb.tscalefactor;                        /* temper. scale */
        vy *= glb.tscalefactor;
        vz *= glb.tscalefactor;
        rx1 = water_rx[i][HY1];                                 /* current r(t) */
        ry1 = water_ry[i][HY1];
        rz1 = water_rz[i][HY1];
        px1 = rx1 + vx * glb.dt;                       /* r(t+dt) */
        py1 = ry1 + vy * glb.dt;
        pz1 = rz1 + vz * glb.dt;

        //aptr = &water[i].atoms[OXY];   /* oxygen */
        om2 = 1.0/(water_mass[i][OXY]);
        vx = water_vx[i][OXY] + water_fx[i][OXY] * om2 * glb.dt;       /* v(t+.5dt) */
        vy = water_vy[i][OXY] + water_fy[i][OXY] * om2 * glb.dt;
        vz = water_vz[i][OXY] + water_fz[i][OXY] * om2 * glb.dt;
        vx *= glb.tscalefactor;                         /* temper. scale */
        vy *= glb.tscalefactor;
        vz *= glb.tscalefactor;
        rx2 = water_rx[i][OXY];                                 /* current r(t) */
        ry2 = water_ry[i][OXY];
        rz2 = water_rz[i][OXY];
        px2 = rx2 + vx * glb.dt;                       /* r(t+dt) */
        py2 = ry2 + vy * glb.dt;
        pz2 = rz2 + vz * glb.dt;

        //aptr = &water[i].atoms[HY2];   /* 2nd hydrogen */
        om3 = 1.0/(water_mass[i][HY2]);
        vx = water_vx[i][HY2] + water_fx[i][HY2] * om3 * glb.dt;       /* v(t+.5dt) */
        vy = water_vy[i][HY2] + water_fy[i][HY2] * om3 * glb.dt;
        vz = water_vz[i][HY2] + water_fz[i][HY2] * om3 * glb.dt;
        vx *= glb.tscalefactor;                              /* temper. scale */
        vy *= glb.tscalefactor;
        vz *= glb.tscalefactor;
        rx3 = water_rx[i][HY2];                                 /* current r(t) */
        ry3 = water_ry[i][HY2];
        rz3 = water_rz[i][HY2];
        px3 = rx3 + vx * glb.dt;                       /* r(t+dt) */
        py3 = ry3 + vy * glb.dt;
        pz3 = rz3 + vz * glb.dt;

        if(DODUMPVERBOSE == 1){
        LOG.lout.printf("i:%d px1:%lg py1:%lg pz1:%lg px2:%lg py2:%lg pz2:%lg px3:%lg py3:%lg pz3:%lg\n",i,px1,py1,pz1,px2,py2,pz2,px3,py3,pz3);
        }
        /* calculate relative vectors (for constraint calculations) */
        rx12 = rx1 - rx2;
        mbound(rx12,glb.box, glb.boxinv);
        ry12 = ry1 - ry2;
        mbound(ry12,glb.box, glb.boxinv);
        rz12 = rz1 - rz2;
        mbound(rz12,glb.box, glb.boxinv);
        rx23 = rx2 - rx3;
        mbound(rx23,glb.box, glb.boxinv);
        ry23 = ry2 - ry3;
        mbound(ry23,glb.box, glb.boxinv);
        rz23 = rz2 - rz3;
        mbound(rz23,glb.box, glb.boxinv);
        rx31 = rx3 - rx1;
        mbound(rx31,glb.box, glb.boxinv);
        ry31 = ry3 - ry1;
        mbound(ry31,glb.box, glb.boxinv);
        rz31 = rz3 - rz1;
        mbound(rz31,glb.box, glb.boxinv);
        px12 = px1 - px2;
        mbound(px12,glb.box, glb.boxinv);
        py12 = py1 - py2;
        mbound(py12,glb.box, glb.boxinv);
        pz12 = pz1 - pz2;
        mbound(pz12,glb.box, glb.boxinv);
        px23 = px2 - px3;
        mbound(px23,glb.box, glb.boxinv);
        py23 = py2 - py3;
        mbound(py23,glb.box, glb.boxinv);
        pz23 = pz2 - pz3;
        mbound(pz23,glb.box, glb.boxinv);
        px31 = px3 - px1;
        mbound(px31,glb.box, glb.boxinv);
        py31 = py3 - py1;
        mbound(py31,glb.box, glb.boxinv);
        pz31 = pz3 - pz1;
        mbound(pz31,glb.box, glb.boxinv);

        /* calculate scalar products (for eqns 3.61a,b,c squared) */
        r12r23 = (rx12 * rx23) + (ry12 * ry23) + (rz12 * rz23);
        r23r31 = (rx23 * rx31) + (ry23 * ry31) + (rz23 * rz31);
        r31r12 = (rx31 * rx12) + (ry31 * ry12) + (rz31 * rz12);
        p12sq  = (px12 * px12) + (py12 * py12) + (pz12 * pz12);
        p23sq  = (px23 * px23) + (py23 * py23) + (pz23 * pz23);
        p31sq  = (px31 * px31) + (py31 * py31) + (pz31 * pz31);
        p12r12 = (px12 * rx12) + (py12 * ry12) + (pz12 * rz12);
        p12r23 = (px12 * rx23) + (py12 * ry23) + (pz12 * rz23);
        p12r31 = (px12 * rx31) + (py12 * ry31) + (pz12 * rz31);
        p23r12 = (px23 * rx12) + (py23 * ry12) + (pz23 * rz12);
        p23r23 = (px23 * rx23) + (py23 * ry23) + (pz23 * rz23);
        p23r31 = (px23 * rx31) + (py23 * ry31) + (pz23 * rz31);
        p31r12 = (px31 * rx12) + (py31 * ry12) + (pz31 * rz12);
        p31r23 = (px31 * rx23) + (py31 * ry23) + (pz31 * rz23);
        p31r31 = (px31 * rx31) + (py31 * ry31) + (pz31 * rz31);
        const1 = r12sq - p12sq;
        const2 = r23sq - p23sq;
        const3 = r31sq - p31sq;

        /* calculate matrix, ignoring quadratic terms (see page 94) */
        mat[1-1][1-1] =  2.0 * ( om1 + om2 ) * p12r12;
        mat[1-1][2-1] = -2.0 *   om2         * p12r23;
        mat[1-1][3-1] = -2.0 *   om1         * p12r31;
        mat[2-1][1-1] = -2.0 *   om2         * p23r12;
        mat[2-1][2-1] =  2.0 * ( om2 + om3 ) * p23r23;
        mat[2-1][3-1] = -2.0 *   om3         * p23r31;
        mat[3-1][1-1] = -2.0 *   om1         * p31r12;
        mat[3-1][2-1] = -2.0 *   om3         * p31r23;
        mat[3-1][3-1] =  2.0 * ( om1 + om3 ) * p31r31;

        /* invert matrix */
        inv[1-1][1-1] = mat[2-1][2-1] * mat[3-1][3-1] - mat[2-1][3-1] * mat[3-1][2-1];
        inv[2-1][1-1] = mat[3-1][1-1] * mat[2-1][3-1] - mat[3-1][3-1] * mat[2-1][1-1];
        inv[3-1][1-1] = mat[2-1][1-1] * mat[3-1][2-1] - mat[2-1][2-1] * mat[3-1][1-1];
        inv[1-1][2-1] = mat[3-1][2-1] * mat[1-1][3-1] - mat[3-1][3-1] * mat[1-1][2-1];
        inv[2-1][2-1] = mat[1-1][1-1] * mat[3-1][3-1] - mat[1-1][3-1] * mat[3-1][1-1];
        inv[3-1][2-1] = mat[3-1][1-1] * mat[1-1][2-1] - mat[3-1][2-1] * mat[1-1][1-1];
        inv[1-1][3-1] = mat[1-1][2-1] * mat[2-1][3-1] - mat[1-1][3-1] * mat[2-1][2-1];
        inv[2-1][3-1] = mat[2-1][1-1] * mat[1-1][3-1] - mat[2-1][3-1] * mat[1-1][1-1];
        inv[3-1][3-1] = mat[1-1][1-1] * mat[2-1][2-1] - mat[1-1][2-1] * mat[2-1][1-1];
        determ = mat[1-1][1-1] * inv[1-1][1-1] + mat[1-1][2-1] * inv[2-1][1-1] + mat[1-1][3-1] * inv[3-1][1-1];
        if(DODUMPVERBOSE == 1){
           LOG.lout.printf("determ:%lg\n",determ);
        }
        if(    Math.abs( determ ) < 1.0e-7){
            if(DODUMPERRORS == 1){
                for (int ii=0; ii < 3; ii++){
                     for (int jj=0; jj < 3; jj++){
                          LOG.lout.printf("mat[%d][%d]=%lg, inv[%d][%d]=%lg\n",ii,jj,mat[ii][jj],ii,jj,inv[ii][jj]);
                     }
                }
             //dumper(-1000.0);
            }
            LOG.lout.printf("zero determinent: %lg\n",determ);
           // finish(holdbuf);
        }
        inv[1-1][1-1] /= determ;
        inv[1-1][2-1] /= determ;
        inv[1-1][3-1] /= determ;
        inv[2-1][1-1] /= determ;
        inv[2-1][2-1] /= determ;
        inv[2-1][3-1] /= determ;
        inv[3-1][1-1] /= determ;
        inv[3-1][2-1] /= determ;
        inv[3-1][3-1] /= determ;

        /* calculate quadratic terms for iterative solution */
        q111 = -r12sq * ( om1 + om2 ) * ( om1 + om2);
        q122 = -r23sq * om2 * om2;
        q133 = -r31sq * om1 * om1;
        q112 =  2.0 * r12r23 * ( om1 + om2 ) * om2;
        q123 = -2.0 * r23r31 * om1 * om2;
        q131 =  2.0 * r31r12 * ( om1 + om2 ) * om1;
        q211 = -r12sq * om2 * om2;
        q222 = -r23sq * ( om2 + om3 ) * ( om2 + om3 );
        q233 = -r31sq * om3 * om3;
        q212 =  2.0 * r12r23 * ( om2 + om3 ) * om2;
        q223 =  2.0 * r23r31 * ( om2 + om3 ) * om3;
        q231 = -2.0 * r31r12 * om2 * om3;
        q311 = -r12sq * om1 * om1;
        q322 = -r23sq * om3 * om3;
        q333 = -r31sq * ( om1 + om3 ) * ( om1 + om3 );
        q312 = -2.0 * r12r23 * om1 * om3;
        q323 =  2.0 * r23r31 * ( om1 + om3 ) * om3;
        q331 =  2.0 * r31r12 * ( om1 + om3 ) * om1;

        l12 = l23 = l31 = 0.0;  /* lagrangian multipliers (* dtsq) */
        kcnt = 0;
        while(true){
            quad1 = q111 * l12 * l12 + q112 * l12 * l23
                  + q122 * l23 * l23 + q123 * l23 * l31
                  + q133 * l31 * l31 + q131 * l31 * l12;
            quad2 = q211 * l12 * l12 + q212 * l12 * l23
                  + q222 * l23 * l23 + q223 * l23 * l31
                  + q233 * l31 * l31 + q231 * l31 * l12;
            quad3 = q311 * l12 * l12 + q312 * l12 * l23
                  + q322 * l23 * l23 + q323 * l23 * l31
                  + q333 * l31 * l31 + q331 * l31 * l12;
            vec1 = const1 + quad1;
            vec2 = const2 + quad2;
            vec3 = const3 + quad3;

            /* new guess for solution of matrix equation */
            l12new = inv[1-1][1-1] * vec1 + inv[1-1][2-1] * vec2 + inv[1-1][3-1] * vec3;
            l23new = inv[2-1][1-1] * vec1 + inv[2-1][2-1] * vec2 + inv[2-1][3-1] * vec3;
            l31new = inv[3-1][1-1] * vec1 + inv[3-1][2-1] * vec2 + inv[3-1][3-1] * vec3;
            if(DODUMPVERBOSE == 1){
               LOG.lout.printf("i:%d kcnt:%d l12new:%24.16lg l23new:%24.16lg l31new:%24.16lg l12:%24.16lg l23:%24.16lg l31:%24.16lg ltolerance:%24.16lg\n",i,kcnt,l12new,l23new,l31new,l12,l23,l31,ltolerance);
            }
            if ((Math.abs(l12new - l12) < ltolerance) && (Math.abs(l23new - l23) < ltolerance) && (Math.abs(l31new - l31) < ltolerance)){
                l12 = l12new;
                l23 = l23new;
                l31 = l31new;
                break;
            }
            kcnt++;
            if (kcnt > MAXKCNT){
            if(DODUMPERRORS == 1){
                LOG.lout.printf("i:%d kcnt:%d l12new:%24.16lg l23new:%24.16lg l31new:%24.16lg l12:%24.16lg l23:%24.16lg l31:%24.16lg ltolerance:%24.16lg\n",i,kcnt,l12new,l23new,l31new,l12,l23,l31,ltolerance);
            }
                //dumper(-2000.0);
                //finish("constraint loop abort");
            }
            l12 = l12new;
            l23 = l23new;
            l31 = l31new;
        }

        /* Calculate constrained positions: eqn 3.60a,b,c */
        px1 += om1 * ( l12 * rx12 - l31 * rx31 );
        py1 += om1 * ( l12 * ry12 - l31 * ry31 );
        pz1 += om1 * ( l12 * rz12 - l31 * rz31 );
        px2 += om2 * ( l23 * rx23 - l12 * rx12 );
        py2 += om2 * ( l23 * ry23 - l12 * ry12 );
        pz2 += om2 * ( l23 * rz23 - l12 * rz12 );
        px3 += om3 * ( l31 * rx31 - l23 * rx23 );
        py3 += om3 * ( l31 * ry31 - l23 * ry23 );
        pz3 += om3 * ( l31 * rz31 - l23 * rz23 );

        /* save results, calculate constrained velocity and kinetic energy */
        //aptr = &water[i].atoms[HY1];   /* 1st hydrogen */
        vx = (px1 - rx1) * kenertmp;    /* constrained velocities [berendsen et al, 1984, step 8c] */
        vy = (py1 - ry1) * kenertmp;
        vz = (pz1 - rz1) * kenertmp;
        hke += 0.5 * water_mass[i][HY1] * (vx*vx + vy*vy + vz*vz);
        water_vx[i][HY1] = vx;
        water_vy[i][HY1] = vy;
        water_vz[i][HY1] = vz;
        water_rx[i][HY1] = px1 * glb.pscalefactor;  /* pressure scaling [berendsen et al, 1984, step 9p] */
        water_ry[i][HY1] = py1 * glb.pscalefactor;
        water_rz[i][HY1] = pz1 * glb.pscalefactor;
        //aptr = &water[i].atoms[OXY];   /* oxygen */
        vx = (px2 - rx2) * kenertmp;    /* constrained velocities */
        vy = (py2 - ry2) * kenertmp;
        vz = (pz2 - rz2) * kenertmp;
        hke += 0.5 * water_mass[i][OXY] * (vx*vx + vy*vy + vz*vz);
        water_vx[i][OXY] = vx;
        water_vy[i][OXY] = vy;
        water_vz[i][OXY] = vz;
        water_rx[i][OXY] = px2 * glb.pscalefactor;
        water_ry[i][OXY] = py2 * glb.pscalefactor;
        water_rz[i][OXY] = pz2 * glb.pscalefactor;
        //aptr = &water[i].atoms[HY2];   /* 2nd hydrogen */
        vx = (px3 - rx3) * kenertmp;    /* constrained velocities */
        vy = (py3 - ry3) * kenertmp;
        vz = (pz3 - rz3) * kenertmp;
        hke += 0.5 * water_mass[i][HY2] * (vx*vx + vy*vy + vz*vz);
        water_vx[i][HY2] = vx;
        water_vy[i][HY2] = vy;
        water_vz[i][HY2] = vz;
        water_rx[i][HY2] = px3 * glb.pscalefactor;
        water_ry[i][HY2] = py3 * glb.pscalefactor;
        water_rz[i][HY2] = pz3 * glb.pscalefactor;

        /* calculate virial constraint contribution */
        glb.wc += l12 * r12sq + l23 * r23sq + l31 * r31sq;      /* [allen-tildsley p.49] */
        /* printf(" gbl.wc:%lg l12:%lg r12sq:%lg l23:%lg r23sq:%lg l31:%lg r31sq:%lg\n",gbl.wc,l12,r12sq,l23,r23sq,l31,r31sq); */
    }
    glb.wc /= (3.0 * glb.dtsq);
    upthermo(glb.ke,hke,glb.cnt);
}
private void evalscalefactors(double hpe, gbl glb){
    double nkbt;
    double sr1,sr3,sr9,tmp1;
    double pelrc,wlrc;           /* pe,virial, long range corrections */

    /* calculate LONG-RANGE CORRECTIONS for LJ [allen-tildsley p.64 f.17] */
    /* these are variables if doing pressure scaling */
    sr1 = glb.sigma / glb.rcut;
    sr3 = sr1 * sr1 * sr1;
    sr9 = sr3 * sr3 * sr3;
    tmp1 = (8.0/9.0) * piconst * (double)NUMMOLECULES * ((double)NUMMOLECULES/glb.boxvolume.cur) * glb.epslon * glb.sigsq * glb.sigma;  /* note density = N/volume */
    pelrc = tmp1 * (sr9 - 3.0 * sr3);
    wlrc = -2.0 * tmp1 * (2.0 * sr9 - 3.0 * sr3);

    /* w[t] wc[t-.5dt],  [allen-tildsley p.98, f.7] add in constraints */
    /* LONG-RANGE CORRECTIONS for LJ [allen-tildsley p.64 f.17] */
    upthermo(glb.pe, hpe + pelrc,glb.cnt);   //upthermo(&gbl.pe, hpe + pelrc);    /* pe[t] */
    upthermo(glb.w, glb.wnc + glb.wc + wlrc, glb.cnt); //upthermo(&gbl.w, gbl.wnc + gbl.wc + wlrc);      /* w[t] */

    upthermo(glb.calctemperature, glb.ke.cur * 2.0/((double)NUMFREEDOM * boltzmanconstant), glb.cnt ); /* temp[t] ke[t-.5dt]   [allen-tildesley p.47 eqn2.50] equipartition theorem with constraints */

    /* calc. pressure from virial theorem, [allen-tildesley p.47 eqn2.54, eqn2.59] */
    nkbt = (double)NUMMOLECULES * boltzmanconstant * glb.calctemperature.cur;
    upthermo(glb.calcpressure, ( nkbt + glb.w.cur)/glb.boxvolume.cur, glb.cnt );

    /* (PV/NKT)  Pressure x V/(NxKxT)  =  1 + W/NKT */
    upthermo(glb.compressibilityfactor, 1.0 + glb.w.cur/nkbt, glb.cnt );

    /* [Berendsen et al 1984. eqn 31, step 2P] */
    tmp1 = 1.0 + (glb.dt/glb.pt) * glb.isothermalbeta * (glb.calcpressure.cur - glb.desiredpressure);
    glb.pscalefactor = Math.pow( tmp1 , (1/3.0));     /* mu */

    /* [Berendsen et al 1984. eqn 34, step 3T] */
    if (glb.calctemperature.cur < 10)
        glb.tscalefactor = 1.0;         /* hack fix: for temperature at t0 */
    else
        glb.tscalefactor = Math.sqrt( 1.0 + (glb.dt/glb.tt) * (glb.temperature / glb.calctemperature.cur - 1.0) ); /* lambda */

   if(FORCE_NVE == 1){
      glb.pscalefactor = 1.0;       /* mu */
      glb.tscalefactor = 1.0;       /* lambda */
   }
}
private double force(gbl glb){
    int i,j,k;
  //  ATOM *aptri;
  //  ATOM *aptrj;
    double rij, rtmp, mei4;
    double qi, qj, rxi, rxj, ryi, ryj, rzi, rzj;
    double rxij,ryij,rzij,rijsq,sr2,sr6,sr12,vij,wij,fij,fxij,fyij,fzij;
    double konst0 = 1.0/(4.0 * piconst * permittivityconst);
    double hpe;

    for (i=0; i < NUMMOLECULES; i++){
        for (j=0; j < NUMATOMS; j++){
            water_fx[i][j] = 0.0;
            water_fy[i][j] = 0.0;
            water_fz[i][j] = 0.0;
        }
    }
    hpe = 0.0;
    glb.wnc = 0.0;

    for (i=0; i < (NUMMOLECULES-1); i++){
        for (j=i+1; j < NUMMOLECULES; j++){
       //     aptri = &water[i].atoms[OXY];   /* oxygen */
       //     aptrj = &water[j].atoms[OXY];   /* oxygen */
            rxij = water_rx[i][OXY] - water_rx[j][OXY];  //aptri->rx - aptrj->rx;
            ryij = water_ry[i][OXY] - water_ry[j][OXY];  //aptri->ry - aptrj->ry;
            rzij = water_rz[i][OXY] - water_rz[j][OXY];  //aptri->rz - aptrj->rz;
            mbound(rxij,glb.box, glb.boxinv);
            mbound(ryij,glb.box, glb.boxinv);
            mbound(rzij,glb.box, glb.boxinv);
            rijsq = rxij * rxij + ryij * ryij + rzij * rzij;
            if (rijsq < glb.rcutsqbig){
                if (rijsq < glb.rcutsq){
                    /* lennard-jones (only for Oxy - Oxy) */
                    sr2 = glb.sigsq / rijsq;
                    sr6 = sr2 * sr2 * sr2;
                    sr12 = sr6 * sr6;
                    vij = sr12 - sr6;
                    hpe += 4.0 * glb.epslon * vij;
                    wij = vij + sr12;
                    wij *= 24.0 * glb.epslon;
                    glb.wnc += wij/3.0;
                    fij = wij / rijsq;
                    fxij = fij * rxij;
                    fyij = fij * ryij;
                    fzij = fij * rzij;
                    water_fx[i][OXY] += fxij; //aptri->fx += fxij;
                    water_fy[i][OXY] += fyij; //aptri->fy += fyij;
                    water_fz[i][OXY] += fzij; //aptri->fz += fzij;
                    water_fx[j][OXY] -= fxij; //aptrj->fx -= fxij;
                    water_fy[j][OXY] -= fyij; //aptrj->fy -= fyij;
                    water_fz[j][OXY] -= fzij; //aptrj->fz -= fzij;
                }

                /* fall thru, O-O rijsq,rxij,ryij,rzij already calculated */
                k = 0;
                qi = qj = 0.0;
                while(true){               /* do the coulomb forces */
                    if (rijsq < glb.rcutsq){
                        if(k == 0){
                            qi = water_charge[i][OXY];
                            qj = water_charge[j][OXY];
                        }else{
                            qi = water_charge[i][nxtitab[k]];
                            qj = water_charge[j][nxtjtab[k]];
                        }
                        rij = Math.sqrt(rijsq);
                        rtmp = rij/ glb.rcut;
                        mei4 = 1.0 + rtmp * (rtmp - 2.0);  /* termination funct */
                        vij = konst0 * qi * qj / rij;  // aptri->charge * aptrj->charge / rij;
                        fij = vij * (1.0/rijsq - 1.0/glb.rcutsq); /* mei4 */
                        glb.wnc += fij * rijsq/3.0;
                        vij *= mei4;        /* prevost et al, brooks et al */
                        hpe += vij;
                        fxij = fij * rxij;
                        fyij = fij * ryij;
                        fzij = fij * rzij;
                        /* printf("i:%d j:%d rij:%lg vij:%lg konst0:%lg fij:%lg fxij:%lg fyij:%lg fzij:%lg rxij:%lg ryij:%lg rzij:%lg\n",
                            i,j,rij,vij,konst0,fij,fxij,fyij,fzij,rxij,ryij,rzij); */
                        water_fx[i][OXY] += fxij; //aptri->fx += fxij;
                        water_fy[i][OXY] += fyij; //aptri->fy += fyij;
                        water_fz[i][OXY] += fzij; //aptri->fz += fzij;
                        water_fx[j][OXY] -= fxij; //aptrj->fx -= fxij;
                        water_fy[j][OXY] -= fyij; //aptrj->fy -= fyij;
                        water_fz[j][OXY] -= fzij; //aptrj->fz -= fzij;
                    }
                    if (k >= 8)
                        break;      /* all done */
                    /* get next combination */
                    //aptri = &water[i].atoms[nxtitab[k]];
                    //aptrj = &water[j].atoms[nxtjtab[k]];
                    k++;
                    rxij = water_rx[i][nxtitab[k]] - water_rx[j][nxtjtab[k]];  //aptri->rx - aptrj->rx;
                    ryij = water_ry[i][nxtitab[k]] - water_ry[j][nxtjtab[k]];  //aptri->ry - aptrj->ry;
                    rzij = water_rz[i][nxtitab[k]] - water_rz[j][nxtjtab[k]];  //aptri->rz - aptrj->rz;
                    mbound(rxij,glb.box, glb.boxinv);
                    mbound(ryij,glb.box, glb.boxinv);
                    mbound(rzij,glb.box, glb.boxinv);
                    rijsq = rxij * rxij + ryij * ryij + rzij * rzij;
                }
            }
        }
    }
    return(hpe);
}
private void savestate(gbl glb){
    int i,j;
    double groo,groh,grhh;
    try{
            FileWriter fout = new FileWriter("Restart_1.rst");
            PrintWriter pout = new PrintWriter(fout, true);
            pout.printf("\n\n %d",1);
    }catch( IOException e ){System.err.println( e );}
}
private void createacfvelbins(gbl glb){
    int i,j,k;
    double vx,vy,vz,tmp;

    /* use oxygen as center of mass (hack) */
    for (j=0; j < ACFNUMREF; j++){
        k = (glb.vacfptr + j * VACFSPACING) % NUMVACFBINS; /* tricky: wrap-around, k goes from 0 to 100 */

        if (k == 0){  /* save velocities */
            for (i=0; i < NUMMOLECULES; i++){
            //    aptri = &water[i].atoms[ OXY ];
                acfvel_vx[j][i] = water_vx[i][3];
                acfvel_vy[j][i] = water_vy[i][3];
                acfvel_vz[j][i] = water_vz[i][3];
            }
            acfvel_validflag[j] = 1;
        }

        if (acfvel_validflag[j] == 1){
            tmp = 0.0;
            for (i=0; i < NUMMOLECULES; i++){
                //aptri = &water[i].atoms[ OXY ];
                vx = water_vx[i][3];
                vy = water_vy[i][3];
                vz = water_vz[i][3];
                tmp +=  (vx * acfvel_vx[j][i]) + (vy * acfvel_vy[j][i]) + (vz * acfvel_vx[j][i]);
            }
            vacf_total[k] += tmp;
            vacf_cnt[k]++;
        }
    }
    glb.vacfptr++;
    if (glb.vacfptr >= NUMVACFBINS)
        glb.vacfptr = 0;
}
private double calcselfdiff(double t,gbl glb){
      int i,j;
  //  ATOM *aptri;
    double ktot,tmp,kdiff;

    /* use oxygen as center of mass (hack) */
    ktot = 0.0;
    for (i=0, j=RXINITWIDTH; i < NUMMOLECULES; i++, j+= (RXINITWIDTH*NUMATOMS)){
      //  aptri = &water[i].atoms[ OXY ];
        tmp = water_vx[i][3] - rxinit[j];
        ktot += tmp * tmp;
        /* printf("i:%d tmp:%lg rx:%lg rx0:%lg ktot:%lg\n",i,tmp,aptri->rx,rxinit[j],ktot); */
        tmp = water_vy[i][3] - rxinit[j+1];
        ktot += tmp * tmp;
        tmp = water_vz[i][3] - rxinit[j+2];
        ktot += tmp * tmp;
    }
    kdiff = ktot/(2.0 * 3.0 * NUMMOLECULES * t); /* Einstein relation, [Allen-Tildesley p.60 eqn2.110] */
    return(kdiff);
}
private void calcpairdist(gbl glb){
    int i,j,jtmp,k;
  //  ATOM *aptri;
 //   ATOM *aptrj;
    double rij, rtmp;
    double rxij,ryij,rzij,rijsq;
    double konst0 = 1.0/(4.0 * piconst * permittivityconst);
    double v, mei4;
    double sr2,sr6,sr12,vij;
    double loc_j2kcal = JOULE2CAL * avogadroconst /(KILO); /* kcal/mole of water molecules */
    int hbondcnt = 0;
//    for (i=0; i < NUMMOLECULES; i++)
//        water[i].kcolor = LOC_BLUE;

    glb.paircnt++;
    for (i=0; i < (NUMMOLECULES-1); i++){
        for (j=i+1; j < NUMMOLECULES; j++){
        //    aptri = &water[i].atoms[OXY];   /* oxygen */
        //    aptrj = &water[j].atoms[OXY];   /* oxygen */
            rxij = water_rx[i][OXY] - water_rx[j][OXY]; //aptri->rx - aptrj->rx;
            ryij = water_ry[i][OXY] - water_ry[j][OXY]; //aptri->ry - aptrj->ry;
            rzij = water_rz[i][OXY] - water_rz[j][OXY]; //aptri->rz - aptrj->rz;
            mbound(rxij,glb.box, glb.boxinv);
            mbound(ryij,glb.box, glb.boxinv);
            mbound(rzij,glb.box, glb.boxinv);
            rijsq = rxij * rxij + ryij * ryij + rzij * rzij;
            if (rijsq < glb.rcutsqbig){
                v = 0.0;
                if (rijsq < glb.rcutsq){
                    /* lennard-jones (only for Oxy - Oxy) */
                    sr2 = glb.sigsq / rijsq;
                    sr6 = sr2 * sr2 * sr2;
                    sr12 = sr6 * sr6;
                    vij = sr12 - sr6;
                    v += 4.0 * glb.epslon * vij;
                }
                /* fall thru, O-O rijsq,rxij,ryij,rzij already calculated */
                k = 0;
                while(true){               /* do the coulomb forces */
                    if (rijsq < glb.rcutsq){
                        rij = Math.sqrt(rijsq);

                        rtmp = rij / glb.bindelta;
                        jtmp = (int)rtmp;
                        if (jtmp >= MAXBIN)
                            jtmp = MAXBIN-1;
                        pairs_bin[jtmp][ binseltab[k] ]++;

                        rtmp = rij/ glb.rcut;
                        mei4 = 1.0 + rtmp * (rtmp - 2.0);  /* termination funct */
                        vij = konst0 * water_charge[i][3] * water_charge[j][3] / rij;
                        v += vij * mei4;   /* prevost et al, brooks et al */
                    }
                    if (k >= 8)
                        break;      /* all done */
                    /* get next combination */
                    //    aptri = &water[i].atoms[nxtitab[k]];
                    //    aptrj = &water[j].atoms[nxtjtab[k]];
                    k++;
                    rxij = water_rx[i][nxtitab[k]] - water_rx[j][nxtitab[k]];  //aptri->rx - aptrj->rx;
                    ryij = water_ry[i][nxtitab[k]] - water_ry[j][nxtitab[k]];  //aptri->ry - aptrj->ry;
                    rzij = water_rz[i][nxtitab[k]] - water_rz[j][nxtitab[k]];  //aptri->rz - aptrj->rz;
                    mbound(rxij,glb.box, glb.boxinv);
                    mbound(ryij,glb.box, glb.boxinv);
                    mbound(rzij,glb.box, glb.boxinv);
                    rijsq = rxij * rxij + ryij * ryij + rzij * rzij;
                }

                v *= loc_j2kcal;
                if ((i == TARGETMOLECULE) || (j == TARGETMOLECULE)){
                    if (v > 0.0)
                        jtmp = LOC_GREEN;
                    else if (v < POTENTIALHYDROBOND){
                        jtmp = LOC_RED;
                        hbondcnt++;
                    }
                    else
                        jtmp = LOC_BLUE;
                    if (i == TARGETMOLECULE)
                        water_kcolor[j] = jtmp;
                    else
                        water_kcolor[i] = jtmp;
                }
                v += EPDISTMIN;
                v *= EPDISTMULT;
                jtmp = (int)v;
                if (jtmp < 0)
                    jtmp = 0;
                if (jtmp >= MAXEPDIST)
                    jtmp = MAXEPDIST-1;
                energypairdist[jtmp]++;
            }
        }
    }

    water_kcolor[TARGETMOLECULE] = LOC_YELLOW;
    upthermo(glb.hbonds, (double)hbondcnt, glb.cnt);
    if (hbondcnt < MAXHBONDBINS)
        glb.hbondbins[hbondcnt]++;
}
private void initlow(int iflag, gbl globStruc){
    int i;

    initthermo(globStruc.w, 0.0,iflag);
    initthermo(globStruc.ke, 0.0,iflag);
    initthermo(globStruc.pe, 0.0,iflag);
    initthermo(globStruc.energy, 0.0,iflag);
    initthermo(globStruc.calcpressure, 0.0,iflag);
    initthermo(globStruc.compressibilityfactor, 0.0,iflag);        /* pv/nkt */
    initthermo(globStruc.calctemperature, 0.0,iflag);
    initthermo(globStruc.hbonds, 0.0,iflag);

    initthermo(globStruc.boxvolume, globStruc.box * globStruc.box * globStruc.box,iflag); /* stillinger   density = 1gm/cm3 */
    globStruc.paircnt = 0;
    globStruc.vacfptr = 0;                    /* for veloc auto correl */
    globStruc.evalcnt = 0;
    globStruc.pscalefactor = 1.0;                /* mu */
    globStruc.tscalefactor = 1.0;                /* lambda */

    for (i=0; i < MAXHBONDBINS; i++)
        globStruc.hbondbins[i] = 0;

    for (i=0; i < MAXBIN; i++){
        pairs_bin[i][0] = 0;
        pairs_bin[i][1] = 0;
        pairs_bin[i][2] = 0;
    }

    for (i=0; i < MAXEPDIST; i++)
        energypairdist[i] = 0;

    for (i=0; i < ACFNUMREF; i++)
        acfvel_validflag[i] = 0;

    for (i=0; i < NUMVACFBINS; i++){
        vacf_total[i] = 0;
        vacf_cnt[i] = 0;
    }
}
private void initthermo(THERMO ptr, double val, int iflag){
    if(iflag == 1)              /* init with val, else no change */
        ptr.cur = val;
    ptr.min = BIGNUMBER;
    ptr.max = -BIGNUMBER;
    ptr.tot = 0.0;
    ptr.rms = 0.0;
    ptr.cnt = 0.0;
}
private void upthermo(THERMO ptr, double val, int cnt){
    ptr.cur = val;
    if(cnt < WAITAVG)
        return;
    if(val < ptr.min)
        ptr.min = val;
    if(val > ptr.max)
        ptr.max = val;
    ptr.tot += val;               /*for average */
    ptr.rms += val * val;         /*for rms fluctuation */
    ptr.cnt += 1.0;               /*for average and rms*/
}
private void initUniverseConstants()   /* SI units */{
    /* FUNDAMENTAL CONSTANTS */
    piconst =  Math.atan(1.0) * 4.0;  /* 3.14159265358979  3.1415926539 */
    boltzmanconstant = 1.380657e-23;    /* joule/kelvin */
    avogadroconst = 6.0221367e23;       /* particles/mole */
    elementarycharge = 1.60217738e-19;        /* coulomb */
    permittivityconst = 8.85418781762e-12;       /* coulombSq/(Newton*MeterSq) */
    deg2rad = piconst/180.0;            /* degrees to radian */
}
public void mbound(double r12, double bx, double ibx){
    double tmp;
    tmp = (r12)*ibx;
    if(tmp < 0.0){tmp -= 0.5;}else{tmp += 0.5;}
    tmp = (double)((int)(tmp));
    r12 -= (tmp*bx);
}
  int RXINITWIDTH = 6;
  double rxinit[] = {
  2.647020211374952e-010,   -4.23539182794927e-010,   6.426252610119555e-010,       -307.8394253804442,        210.0946626679413,        1280.499133389315, /* 0 HY1 */
  3.570092548791954e-010,   -4.03150064824816e-010,   6.752390382748862e-010,         277.944046560782,       -177.0323434563364,       -116.6891508061831, /* 0 OXY */
   3.70615558442565e-010,  -3.040957587463092e-010,   6.770033938646047e-010,        162.9688916819176,       -180.1563640372302,        1188.683246382176, /* 0 HY2 */
  3.104596672017118e-009,  -1.306136279464044e-010,  -1.530129316886179e-010,       -1057.111920344608,       -1484.047866082394,        470.9687629860098, /* 1 HY1 */
  3.060686205442587e-009,  -1.285550219289682e-010,  -2.428329353986648e-010,        204.6702639661331,       -295.5558199194026,       -128.0030589404391, /* 1 OXY */
   3.06069679734584e-009,  -2.205419650016775e-010,  -2.820554058964698e-010,        2404.248660039809,        547.8403997327637,       -2167.271907730412, /* 1 HY2 */
 -7.481646719775176e-010,   3.253879033644156e-010,  -6.018661212062809e-010,       -504.3636080184531,       -157.7756790371572,        986.4002499277662, /* 2 HY1 */
 -8.318767059175719e-010,    3.58640571587405e-010,  -6.453006170554563e-010,        190.8787781458412,        428.7448762610808,        86.04950377117511, /* 2 OXY */
 -8.445744954751062e-010,   4.554029737554326e-010,  -6.234876927437168e-010,       -1602.877933702088,       -537.6729816639646,        3510.922159566834, /* 2 HY2 */
  1.798129369372563e-009,  -2.948553679951742e-011,  -1.346511717678368e-010,       -705.4994359979662,       -197.6673592557205,       -319.7396183166191, /* 3 HY1 */
  1.808846540628275e-009,    2.28888879542027e-011,  -2.191618922238448e-010,       -72.25813433395456,       -178.3423787716366,       -228.6709965233832, /* 3 OXY */
  1.732694231328247e-009,   2.721234986652483e-012,  -2.807580898504524e-010,        819.2334468960291,       -789.5877464113503,       -1138.807367398471, /* 3 HY2 */
  6.578367199410607e-010,  -7.823778183594556e-010,  -3.162504017346842e-009,          1085.4008312179,        1106.272678924233,        1573.683736468174, /* 4 HY1 */
  7.400999093376439e-010,  -7.693929923716891e-010,  -3.217858890342546e-009,        268.5975087160514,       -169.5967887796934,        39.54715098833799, /* 4 OXY */
  7.979586290082948e-010,  -6.995663345971595e-010,  -3.175709081065222e-009,        74.03091403182979,        2109.320692927934,       -3368.768003586212, /* 4 HY2 */
 -3.149331577721648e-010,  -3.804738098299126e-010,  -1.349750334924922e-010,       -457.6827390521798,        185.3165813919961,       -1214.159869395638, /* 5 HY1 */
 -2.476649615292756e-010,  -3.325113156562643e-010,  -1.913185240381832e-010,        198.5557842672411,        353.1246867825471,       -293.6288539055841, /* 5 OXY */
 -1.819583146153561e-010,  -2.857305088830524e-010,  -1.322068274885032e-010,         168.161103861188,       -734.3095614180086,        609.2012441808306, /* 5 HY2 */
 -7.156178194434096e-011,   -1.18425525083392e-009,   6.477017938322807e-011,        1723.188503340178,       -2612.063249581858,        817.1621783405859, /* 6 HY1 */
 -7.649476024783038e-011,  -1.128614213489374e-009,   1.477143816416273e-010,        328.2969904354976,       -64.74760032377489,       -939.8723502773618, /* 6 OXY */
  1.377677154898789e-011,  -1.124556637955562e-009,   1.905464563969512e-010,       -361.0354745430602,        1315.708324568281,        406.6485160804452, /* 6 HY2 */
  1.268345603229967e-009,   1.338799189050857e-009,  -1.539855432054624e-009,        817.3025698003245,       -861.1430707720009,        660.2521040302057, /* 7 HY1 */
  1.245288823137777e-009,   1.249453604543757e-009,  -1.578400878230088e-009,        309.8176621319864,        -232.927620051794,       -505.3858648442639, /* 7 OXY */
  1.201070284877198e-009,   1.193363552121108e-009,   -1.50841068682213e-009,        -1668.68569775403,         -254.18787226412,        -1752.86399167031, /* 7 HY2 */
 -6.130089841871545e-011,   1.470218774358181e-009,   1.573658845338955e-009,         296.842739661919,        1423.413346569967,        1583.280569180126, /* 8 HY1 */
  2.925645012144076e-011,   1.434124566472802e-009,   1.595942803302132e-009,        265.4742233954022,        263.9023394638062,       -119.7188505003856, /* 8 OXY */
  2.761716915644702e-011,   1.394284423708232e-009,   1.687649268627205e-009,        1223.239704708521,       -626.3658752981908,        -484.334983615991, /* 8 HY2 */
  5.855817125997997e-010,   6.652514811852557e-010,  -2.090674561114224e-009,        2626.155432913405,        3107.149493221515,       -198.1302457444028, /* 9 HY1 */
  5.188326308131087e-010,   7.289800726579459e-010,  -2.129187243564823e-009,       -260.0818985282968,        269.8389937791949,        2.599544731373004, /* 9 OXY */
  4.636143909764827e-010,   7.679786410545183e-010,  -2.055498318668589e-009,        1760.097919081826,          1762.7689948632,        749.5996524425905, /* 9 HY2 */
  1.925321876846785e-009,   4.557614534240676e-010,   8.154373475801251e-010,         1112.52862677359,        363.2391965084398,        1290.911264973627, /* 10 HY1 */
  1.903629352727708e-009,   3.689486431225354e-010,   8.600800442380907e-010,       -355.2426707321878,        251.6661814544666,         377.539714008253, /* 10 OXY */
  1.988023431233024e-009,   3.259818236616834e-010,   8.921968924436324e-010,       -1078.065333215776,       -1499.953330869414,       -37.16198672042686, /* 10 HY2 */
  3.179727496813718e-010,    2.13265975067208e-009,   7.292606403233552e-010,       -5.394154636037759,        2276.729321865576,        1680.205594141626, /* 11 HY1 */
  2.861299902824541e-010,   2.181493542310982e-009,   6.480122589977442e-010,        178.1588235222796,        203.7017057560245,        343.4634473410956, /* 11 OXY */
  2.058616544128431e-010,   2.236349512629832e-009,   6.714177450448441e-010,        1299.335613099777,        2346.108264467477,       -757.3098563240175, /* 11 HY2 */
   1.00820007320076e-009,   9.328700491407128e-010,   -1.68587111417644e-009,       -374.2283897780566,        299.1193321859664,       -38.08617063632348, /* 12 HY1 */
  1.072049016752039e-009,   9.533641015142932e-010,  -1.760055386260862e-009,       -545.1499403520762,        457.7922749920419,       -141.5795953237732, /* 12 OXY */
  1.051674333271123e-009,   1.043807854148097e-009,  -1.797536046374705e-009,       -621.1905503762742,        545.0270841878763,        109.7502144062855, /* 12 HY2 */
  1.311088774927068e-009,    1.12255948317359e-009,   1.754023821135472e-009,       -3226.476733832097,        -467.672927125692,        496.6772241339049, /* 13 HY1 */
  1.318833659806784e-009,   1.046819523394067e-009,   1.689189607553461e-009,       -181.2360414577324,       -255.1552946961042,        576.2296020541294, /* 13 OXY */
  1.388967790383473e-009,   1.068194340660211e-009,   1.621187247924605e-009,         379.032739505529,        1614.486800351965,        1722.897816749794, /* 13 HY2 */
  2.589464701262465e-009,  -4.812597702675014e-011,  -7.054155670272409e-010,        564.7083242937288,        -226.838107722147,        231.0277247737014, /* 14 HY1 */
   2.49315101837446e-009,  -4.834809944699665e-011,  -6.785152854268679e-010,        531.7858333978959,        170.7143978523514,        118.0323996188089, /* 14 OXY */
  2.466017220867907e-009,    4.38824469139685e-011,  -6.509965036436694e-010,        1221.806900095697,        25.38984039953347,        1302.730358925814, /* 14 HY2 */
 -3.610144153657254e-010,   3.637204715645294e-010,  -2.368777914005818e-009,        52.85844431308456,       -906.4768677504172,       -466.3224787280224, /* 15 HY1 */
 -2.903990184110538e-010,   3.228102309345609e-010,  -2.310986756165421e-009,        74.06967692032939,        -172.920449873747,        30.43970659815221, /* 15 OXY */
 -2.048960647137508e-010,    3.16254564315326e-010,  -2.362428647112543e-009,       -205.3515287934076,       -1730.724919379334,       -248.0198146314548, /* 15 HY2 */
  2.635412851209044e-009,   4.997582460032029e-010,  -2.079855574033384e-009,        502.8858667844463,       -624.0849307431072,        1272.194409919433, /* 16 HY1 */
  2.706489508255334e-009,   5.665579300948221e-010,  -2.057812184357918e-009,       -461.4573290880925,        765.4183629022994,        216.0031215382936, /* 16 OXY */
  2.670007901808499e-009,   6.341375721769695e-010,  -1.993764667879072e-009,       -1209.914784442494,        805.0609936321378,         -249.11159036954, /* 16 HY2 */
  1.104110298731588e-009,   2.459668166227647e-010,  -1.002034926727468e-009,        295.3177311816426,        495.4231253477956,        109.6392306938665, /* 17 HY1 */
  1.078772842324828e-009,    3.11261200296588e-010,   -1.07341179107084e-009,         433.738691032052,       -67.37218776517923,        -456.633179198268, /* 17 OXY */
   1.15457765882514e-009,     3.2127497225657e-010,  -1.137858312827164e-009,        1495.562418896261,        1272.157698850064,        981.1149679321044, /* 17 HY2 */
 -6.043153125035568e-010,   1.050734798868464e-009,  -1.322404490058517e-009,        675.2434737722112,        429.4981550803638,        1771.709276961894, /* 18 HY1 */
 -5.604921672834106e-010,   9.608622536738474e-010,  -1.323972151774546e-009,       -169.8323904159797,        35.42992873264437,        214.0229057101933, /* 18 OXY */
 -5.937172249350247e-010,   9.062327465767778e-010,  -1.247084490226165e-009,        1065.004531145403,       -141.3350286471393,        627.6488816500565, /* 18 HY2 */
  6.388120943853174e-010,   1.851135533315215e-009,   2.106162933781083e-009,         148.462384789013,       -165.8513921779868,        777.6262430185124, /* 19 HY1 */
  5.537527437483947e-010,   1.831087840034008e-009,   2.154773599928877e-009,       -447.1980141228492,        484.9933811839689,        10.77943282587353, /* 19 OXY */
  4.768266288346546e-010,   1.871788913777376e-009,   2.105521235108323e-009,        328.3040880676358,         1249.41558600468,       -576.5269392179716, /* 19 HY2 */
 -2.237831562478618e-010,   6.656256993003889e-010,   2.155158754217182e-010,        613.5438679348141,       -375.1335932779279,        91.67530785445874, /* 20 HY1 */
 -2.511952592383288e-010,   6.183793394763206e-010,   2.992795793334666e-010,       -385.6640089964318,        138.2689245474423,        58.03058604217195, /* 20 OXY */
 -3.485951066628575e-010,   5.961948649351541e-010,   2.946840585783181e-010,       -369.2880908650895,        254.0687354114586,       -898.5900737390076, /* 20 HY2 */
  1.067565220485416e-009,   3.158271648691519e-009,   -5.56523919593187e-010,       -622.2574132760269,        1405.882516814626,        189.6788778138982, /* 21 HY1 */
   1.08409542222768e-009,    3.15834832976023e-009,  -4.578996500122431e-010,       -568.0860399754914,       -64.16545384788579,        187.2277216695993, /* 21 OXY */
  1.102072628271756e-009,   3.064942006093939e-009,  -4.270438764520114e-010,       -2260.550570231968,       -769.6656543803545,       -925.1256039605608, /* 21 HY2 */
  1.237938711853399e-009,   6.928908622381764e-010,  -6.947300129231321e-010,        1338.361092677932,        438.0437701709516,       -408.0037082525318, /* 22 HY1 */
  1.238118991632351e-009,   7.619242019704147e-010,  -6.223810979028122e-010,        350.3779878181652,        351.9474515530078,       -319.9657846205199, /* 22 OXY */
  1.197281355493286e-009,   8.464448658162105e-010,  -6.568563768377392e-010,       -1225.236272334038,       -307.1573797732018,       -90.98920284067932, /* 22 HY2 */
 -6.936734405024524e-010,  -6.457319631499225e-010,  -2.325428165573162e-013,        740.3279209025833,       -2133.794957619114,        257.2771616305026, /* 23 HY1 */
 -5.985387107800908e-010,  -6.228115858770441e-010,    2.03596767414379e-011,        337.9581950505774,       -59.75359399473103,        -136.263206599454, /* 23 OXY */
 -5.833534127679649e-010,  -6.290422883421827e-010,   1.190034047279614e-010,        997.0231816603302,       -184.5385036084064,       -244.4321341221825, /* 23 HY2 */
 -1.462575818265162e-009,   1.639294867366234e-009,  -5.911468666015967e-011,       -105.9304880822402,       -212.5215093601852,       -774.7600798643872, /* 24 HY1 */
 -1.373305547244683e-009,   1.645070619723688e-009,  -1.442174094570054e-011,       -306.2973454479928,        585.6478881956475,       -473.3974185911582, /* 24 OXY */
 -1.385491326923056e-009,   1.672598277585554e-009,   8.093933510878904e-011,       -796.8481111021204,        1474.645704749928,       -789.7425635700655, /* 24 HY2 */
  -1.84860218099666e-009,  -2.150012490946162e-010,  -1.900160294032016e-009,        304.9572653634956,        -246.195292426008,        176.3629049042604, /* 25 HY1 */
 -1.943994063386303e-009,   -2.37818356445314e-010,  -1.880672647013752e-009,        181.1628066053009,        63.28635739630303,       -65.08033957044614, /* 25 OXY */
 -1.990610909653497e-009,  -1.574394314664282e-010,  -1.843711883177743e-009,        144.1636580278817,        435.9116379704133,       -916.2485344502736, /* 25 HY2 */
  1.623705190150359e-009,   2.932819121624918e-009,    6.91266325498596e-010,         296.318774697601,        823.0838257356372,        2174.732846354112, /* 26 HY1 */
  1.678592599054272e-009,   3.013973494163364e-009,   6.712328494719911e-010,        87.62815382400329,        471.6217570758388,        124.6834167470253, /* 26 OXY */
  1.742502113453451e-009,   2.993518238838952e-009,   5.970900379601668e-010,       -1044.009720910197,       -710.7022936476646,       -535.0721100228061, /* 26 HY2 */
  2.241891918571164e-009,   2.822989055254446e-009,  -3.267048332710515e-010,        235.2949789679203,        1114.543578661851,        48.15474382616166, /* 27 HY1 */
  2.268035816996722e-009,   2.857294166435478e-009,  -2.364847915842254e-010,        381.5614574366876,        275.3655351782608,         327.083609316739, /* 27 OXY */
  2.366869613271903e-009,   2.847807657386461e-009,  -2.245731644955608e-010,        431.9380411657154,        372.6858751280722,       -10.74963630999923, /* 27 HY2 */
  1.889673710593929e-009,  -6.452839493134967e-010,   8.967334989288953e-010,        -1811.52241890417,        448.6277354415244,        394.1881740328026, /* 28 HY1 */
  1.876505439644854e-009,  -6.016670316456018e-010,   8.077157415241814e-010,       -377.0907808241863,       -28.30175003161622,       -58.68456888683651, /* 28 OXY */
  1.799610205220507e-009,  -6.447892159978941e-010,   7.605182824875455e-010,        -1352.98599192207,        1527.043515437306,         92.2397346948438, /* 28 HY2 */
  1.502192593830125e-009,   1.971098500174092e-009,   5.693545999136475e-010,       -1156.411683991618,        2398.604486516262,        144.0063752857026, /* 29 HY1 */
  1.414139617884367e-009,   1.926570534221792e-009,   5.856006355915385e-010,        25.41447343625599,        -209.965996184046,       -468.2642449178191, /* 29 OXY */
  1.348455328805354e-009,   1.956326906980995e-009,   5.163176474313828e-010,        391.3648106600292,       -4429.372510013849,       -2710.263219075101, /* 29 HY2 */
 -2.139055428267749e-009,  -1.781366370791936e-010,  -3.130275597060183e-010,        588.1763435383206,        464.5955839855968,         1690.20885172341, /* 30 HY1 */
 -2.109615556275211e-009,  -9.034699120418984e-011,  -3.507935763639203e-010,        260.4439199481796,       -372.3745775373306,       -549.4150993884831, /* 30 OXY */
 -2.188414228462727e-009,  -4.276246508490854e-011,  -3.898637969351986e-010,       -522.5339560018408,       -218.9210881102362,        1193.123508671562, /* 30 HY2 */
  5.863353762272544e-010,   1.060809616477346e-011,  -1.382362332886561e-009,       -1318.545900786021,       -1564.401456199539,       -211.7409583183514, /* 31 HY1 */
  5.955356344568872e-010,   1.437133527382213e-011,  -1.282857592899377e-009,       -72.27073469474828,        344.5102315194989,       -386.0327926379848, /* 31 OXY */
  6.906556059985751e-010,   3.380734960422956e-011,  -1.258890274299706e-009,         529.020337758291,       -645.3405457657585,       -1930.807788152028, /* 31 HY2 */
  2.391248737356552e-010,    2.45601413766607e-009,   1.391364832911911e-009,        144.4836205804627,        610.4349099179996,       -353.0056070301634, /* 32 HY1 */
  1.724462172728764e-010,   2.461116856988606e-009,   1.317014856262055e-009,       -261.9768333803936,       -169.7324920146134,       -44.94809237164086, /* 32 OXY */
  2.186404209571344e-010,    2.48841775995551e-009,   1.232630282771582e-009,       -717.6786716434015,       -1077.662891109646,       -592.0972316239196, /* 32 HY2 */
  1.064862458104628e-009,  -6.842297197044991e-011,  -3.687511968460199e-009,        98.62650464236736,       -448.0119878713336,        660.5858241819753, /* 33 HY1 */
  1.027279356974012e-009,   2.099983222808268e-011,  -3.663199983788498e-009,        266.4484825559284,       -223.5708401432755,         98.5469699199386, /* 33 OXY */
  1.081269803005565e-009,   9.287499945566396e-011,  -3.707006285542033e-009,         431.839491117758,       -598.0015811879351,       -313.8857259785902, /* 33 HY2 */
  3.468370145509016e-011,   6.985570835505859e-010,   9.126384146866706e-010,       -1514.446564995517,        463.9597433273834,        408.7161297204383, /* 34 HY1 */
  6.577476788236907e-011,   6.809442578882495e-010,   1.006036133934001e-009,       -820.4908580325771,        -414.297788374152,        15.85207665604241, /* 34 OXY */
  6.973823349967092e-011,   7.670977750361837e-010,    1.05665097783717e-009,        173.5246310392643,       -804.3299433934376,        609.2744588082306, /* 34 HY2 */
   2.57850273724934e-009,   2.325237590473092e-009,   3.728624670211366e-010,        661.1977247866261,        680.8055949680705,        1940.167338148741, /* 35 HY1 */
  2.499087869633524e-009,   2.385051463351437e-009,   3.836132536368106e-010,        224.6765590039129,        417.8641929213098,        250.9477737239974, /* 35 OXY */
  2.416452098469945e-009,   2.329452533067676e-009,    3.92562458936849e-010,        752.6878500477918,        19.15954592597214,        2850.500007758126, /* 35 HY2 */
  7.150835698732137e-010,  -2.476924606071655e-010,  -2.274276787566994e-010,       -737.3427494129874,         1089.04182058783,        536.1151004143935, /* 36 HY1 */
  7.156745006715277e-010,   -3.09437545398306e-010,  -1.487689519367929e-010,       -582.0065786381651,       -234.7951786675656,       -495.2027383433866, /* 36 OXY */
  6.446229391531769e-010,  -2.817089898868176e-010,  -8.409436020980317e-011,       -1714.232990778636,       -2461.261649145684,       -760.1049094759984, /* 36 HY2 */
  2.350116563760696e-010,   -1.33968576504459e-009,  -8.633285129215004e-010,        -155.927263759835,        290.7205703874662,        529.9366918228052, /* 37 HY1 */
  2.831955593563713e-010,  -1.427056540634416e-009,  -8.700112568502662e-010,       -560.4568473793815,        68.10061923008972,        515.7677272717293, /* 37 OXY */
  2.585590119842931e-010,  -1.484749062325459e-009,  -7.921356286055313e-010,       -1277.736141431347,       -56.60838365774602,        198.4876542056909, /* 37 HY2 */
 -1.041902725988141e-010,  -1.156973151916991e-009,   4.664284762276847e-010,       -115.7291629023958,         644.171344830111,       -1367.003034059859, /* 38 HY1 */
 -6.684270556248244e-012,  -1.169793601873937e-009,    4.84545198574094e-010,       -513.2652374236667,        -365.047773373869,        104.5096673766757, /* 38 OXY */
  4.636842686566115e-011,  -1.136361783700956e-009,   4.066495964523526e-010,         1045.16200635312,       -149.1490105678921,        1246.446094323818, /* 38 HY2 */
  1.604327141692017e-009,  -6.525275174831216e-010,   9.303763923527264e-010,        671.2297727259827,         16.2632242872167,       -363.4844289363232, /* 39 HY1 */
  1.587603876292642e-009,  -5.539913849768536e-010,   9.270653524721098e-010,        134.1035053798252,       -62.04862219826831,       -12.67830621144798, /* 39 OXY */
  1.659254654946046e-009,  -5.098187179981207e-010,   8.730751551178862e-010,        2000.132998852066,        982.4072348130846,        3247.873381249426, /* 39 HY2 */
  3.555517466192431e-009,  -1.510270776675502e-009,  -1.357416359015376e-010,        99.43052264225576,        1337.178602478115,        594.3683718383996, /* 40 HY1 */
   3.60360757443565e-009,  -1.484345250014712e-009,  -2.194984777505313e-010,       -49.63653955966336,        256.7144701516345,        170.2522447718176, /* 40 OXY */
  3.539890858896216e-009,   -1.48648910282875e-009,  -2.965412194980798e-010,       -335.2372539798932,       -128.0982223800674,        416.2194016735856, /* 40 HY2 */
  4.948855729195795e-010,  -1.527938084220128e-010,  -7.650327734162194e-010,       -737.7277367657667,        325.1902138856057,        444.4199130425276, /* 41 HY1 */
   4.32771368896205e-010,  -2.168207754803658e-010,  -8.102253961311279e-010,       -204.9402529876506,        218.7350600872013,       -140.5654011518836, /* 41 OXY */
  4.264230219381144e-010,  -3.007894830153544e-010,  -7.562906703433394e-010,        383.4949350236292,       -274.5097526230253,       -834.2530190562037, /* 41 HY2 */
 -8.949450859193647e-010,   9.918269765370892e-010,   1.142646980247503e-009,       -293.4637064806361,       -465.8179903333627,       -496.6494934944761, /* 42 HY1 */
 -8.369408291967961e-010,   9.104619609787602e-010,   1.138743078960926e-009,        574.1057474465171,         120.302014764943,        85.88759208823794, /* 42 OXY */
 -8.081490240141211e-010,   8.945531473391015e-010,   1.044308225210302e-009,        1988.908237755619,         618.257179170024,        427.0892721557454, /* 42 HY2 */
  2.277085185439922e-009,   1.369957282618736e-009,   1.989456117456022e-009,       -1811.037016999869,        579.5267691084078,        1613.772197236138, /* 43 HY1 */
  2.219006543613344e-009,   1.328378620846942e-009,   1.919469788442082e-009,       -11.12123551409118,        70.77180741393494,        404.6295936279468, /* 43 OXY */
  2.201458704271276e-009,   1.232826550532398e-009,   1.943173695146053e-009,        2417.078237120017,       -821.9726537119243,       -1295.431478136699, /* 43 HY2 */
 -3.949991120021804e-010,  -9.320914539761986e-011,   1.278741681494024e-009,        960.7195225576232,         157.829905373721,        133.0649998577078, /* 44 HY1 */
  -3.63976448430147e-010,  -8.786537725446408e-011,   1.183825726634315e-009,       -621.9071454479645,        315.5533356268698,       -382.6897435199766, /* 44 OXY */
 -3.191031459170511e-010,   1.476442557846164e-013,   1.168330734985176e-009,         362.851839210517,       -221.2068984759796,       -600.7508693759322, /* 44 HY2 */
  1.680482616969971e-009,   1.882291563211586e-009,    1.56759999679122e-010,        2125.104596510275,       -1935.303475340455,        36.94243310394886, /* 45 HY1 */
  1.682860020996493e-009,   1.878699504936054e-009,   5.685281741758353e-011,        124.0064620042514,        558.1633497372898,       -125.9705260425028, /* 45 OXY */
  1.610787226481694e-009,    1.93720670334566e-009,    1.96719080004264e-011,       -162.3692590398742,        2083.468124364245,        2757.251677643801, /* 45 HY2 */
  2.036216730497625e-009,  -8.515386629476669e-010,    1.87276391668008e-009,        1991.378867469651,       -1551.289493381258,        757.4954753385254, /* 46 HY1 */
  2.042369959749072e-009,  -8.118061793358602e-010,   1.964325194062048e-009,        210.6895753621952,        358.0831956167871,        68.51082080238973, /* 46 OXY */
  1.951595516929226e-009,  -8.104880318642942e-010,   2.006256841036516e-009,        32.15970193968841,       -2513.554353603268,       -177.9849629819224, /* 46 HY2 */
  1.606510173157042e-009,   2.109578715174279e-009,    2.98052171477186e-009,        1916.253449851843,         769.626178521517,        803.6007338479081, /* 47 HY1 */
  1.628175750049023e-009,    2.01338524931412e-009,   2.963865787786302e-009,       -22.91242105436194,        382.8758927520873,        454.2783318723936, /* 47 OXY */
  1.726975598496478e-009,   1.999735625188378e-009,   2.971096120281167e-009,        41.86300659798364,       -1029.592485553642,       -2687.028974807158, /* 47 HY2 */
 -5.284778566845847e-010,  -6.653684500993083e-011,   9.162601058161558e-010,       -989.5010119114048,        27.22332366717076,        871.4411621459767, /* 48 HY1 */
 -5.376256674095656e-010,   -1.55590140224485e-010,   9.608227378116986e-010,        501.6852076416771,        -187.202358422376,        761.8466960541825, /* 48 OXY */
 -4.608774075023545e-010,  -1.695982852680417e-010,   1.023380522888016e-009,         1753.49331300358,        197.3174609325722,       -672.7382412419263, /* 48 HY2 */
  2.618292428232536e-009,  -2.034446972011998e-010,   2.693472927969532e-009,       -781.9805671313666,        363.6584479423587,        815.0304568415032, /* 49 HY1 */
   2.69798384166645e-009,  -1.476656243650005e-010,   2.716667187139424e-009,       -460.1635033804915,       -11.52691036540781,        614.6599444305623, /* 49 OXY */
  2.714543515494572e-009,  -8.125284716222279e-011,    2.64376231606128e-009,       -322.1971724552672,       -98.59219857974376,        566.5863159663062, /* 49 HY2 */
  1.545182415740212e-009,   7.120674419400087e-010,   1.126977161992722e-009,        2445.568769413784,        2215.050635754364,        86.72465522816043, /* 50 HY1 */
  1.551795731148587e-009,   6.961423735229996e-010,   1.028475099805052e-009,        311.9375464320975,        347.0668172175356,         225.017964741803, /* 50 OXY */
  1.579391385463025e-009,   6.014840575389742e-010,    1.01179328101851e-009,       -747.9414832878123,       -175.4104342101801,        1395.002151264412, /* 50 HY2 */
 -1.431340779434842e-009,  -5.220957531864366e-010,   8.208026142243435e-010,        -840.742400827784,       -383.2957169545034,       -307.1766664458179, /* 51 HY1 */
 -1.432934882480279e-009,  -6.030287569120042e-010,   8.795163923640494e-010,       -603.8445445120066,       -377.7039528991167,       -292.7969494233198, /* 51 OXY */
 -1.490235374888421e-009,  -6.734712645993289e-010,   8.376295605370117e-010,       -470.0605931178205,        -519.959502425986,       -236.8214975175351, /* 51 HY2 */
  9.143063546303718e-010,  -1.145713675737537e-009,  -2.794255695949256e-010,        -2329.55092831729,        1394.653432811914,         283.388624313264, /* 52 HY1 */
  9.506605659499386e-010,  -1.075572024826614e-009,  -3.407323402391558e-010,        437.3779246823867,       -61.07451508400766,        218.7554782888184, /* 52 OXY */
  8.754235742773391e-010,  -1.020072656569349e-009,  -3.762177626983103e-010,        2430.519515054407,        1781.835352611765,       -1190.740093282668, /* 52 HY2 */
   1.73642213996297e-009,   1.298594242742434e-009,  -5.221738883806631e-010,        362.6222540703482,          540.46032266006,         756.438259258545, /* 53 HY1 */
  1.766727931522874e-009,   1.238848562568379e-009,  -4.479309580878034e-010,        117.6216498840499,        107.0172216639738,        508.6827447926304, /* 53 OXY */
  1.838796111795058e-009,   1.284476125311256e-009,  -3.957361028777212e-010,        772.8475633077576,       -674.4324124559656,        292.3089957126338, /* 53 HY2 */
   5.29708816625912e-010,    -3.0720971716807e-010,   8.411990346147505e-010,       -566.0119439699157,       -340.4208663504376,       -1065.596338403459, /* 54 HY1 */
  6.250807879209373e-010,  -3.288878681090164e-010,   8.203604918854127e-010,       -140.4934466366944,         206.188316086432,        285.5862478034894, /* 54 OXY */
  6.289792405546944e-010,   -3.88165421335193e-010,   7.399181289751884e-010,         1362.14701495058,        558.7876548767197,        91.05842972242188, /* 54 HY2 */
 -6.515512139706017e-010,  -1.094998807995077e-009,   7.462397485994014e-010,        730.5735305120223,       -701.3501991100494,        566.2958247386376, /* 55 HY1 */
 -6.606061973782259e-010,  -1.018911662899132e-009,   8.104951267917163e-010,        112.4483133585353,       -502.2758154640641,        245.4983976211691, /* 55 OXY */
 -7.570558533442325e-010,  -1.006073347042089e-009,   8.335741542841231e-010,        60.99348762896395,       -1554.343330178024,        629.3177573968273, /* 55 HY2 */
 -1.459732892256662e-009,  -7.726956017651909e-010,   6.492305549009582e-010,         759.770449252042,       -447.2244969224805,        30.96709115902066, /* 56 HY1 */
 -1.365004347562644e-009,  -7.915146445886882e-010,   6.233009107197134e-010,        739.0246179743257,       -456.7343315688402,       -37.97233488787048, /* 56 OXY */
 -1.353937643440504e-009,  -8.894263439462277e-010,   6.062472545488755e-010,        874.1953389959768,       -553.5801988695084,        599.4135693813673, /* 56 HY2 */
  1.772319362354509e-009,  -2.831356037205673e-009,   8.917803667351028e-010,       -1213.418748518098,       -1509.359003842488,         548.133109234514, /* 57 HY1 */
  1.787063144598009e-009,  -2.863789013855176e-009,   7.983420292562639e-010,        135.7914216943683,       -475.3013614666343,        394.3051046543564, /* 57 OXY */
   1.69959444651377e-009,  -2.866989060597636e-009,   7.499789590859895e-010,         826.259169299913,        31.98915647556248,        -900.491000820276, /* 57 HY2 */
  2.123315405223846e-009,  -6.249484644429659e-010,  -1.480976965205404e-009,        300.9584826379013,        221.3360627764129,        31.84553449009166, /* 58 HY1 */
  2.048349344478992e-009,  -6.294609676908517e-010,  -1.414948736008814e-009,        38.81319332062696,        456.2677313624025,       -248.9613054037032, /* 58 OXY */
  2.049978279153726e-009,  -5.479450639344114e-010,  -1.357047924213152e-009,        -382.826083703407,        813.7699039657397,       -738.0578869640041, /* 58 HY2 */
  4.996941490396659e-010,  -1.863706296002302e-010,   3.208056745987277e-010,        380.7963288680092,        530.0598317850217,       -258.7155843425334, /* 59 HY1 */
  4.594146110105037e-010,  -2.756372586797432e-010,    3.41030109051585e-010,        457.1717744529855,        594.6977503814536,         181.210853375051, /* 59 OXY */
   3.64337745377993e-010,  -2.636046987993254e-010,   3.695889452362354e-010,        227.0445630303407,        596.1603688600329,       -579.9991803824156, /* 59 HY2 */
  5.343193001475629e-010,   2.887961438883259e-009,   8.223960961899113e-010,       -1979.208018049014,       -310.8024785224408,        37.04722927108706, /* 60 HY1 */
  5.304894131651865e-010,    2.80061310390422e-009,   8.709306235090533e-010,         204.107111100834,       -398.0994750400946,         76.8255069490999, /* 60 OXY */
  5.754928489869172e-010,    2.80983508498927e-009,   9.597543018223734e-010,        1376.516883743512,        324.3844577196707,         -585.62576606866, /* 60 HY2 */
  1.658696578638214e-009,   3.399374284402416e-009,   2.580441846305412e-009,        685.2812751059344,        560.9328337516747,        160.2165263355387, /* 61 HY1 */
  1.611266142338734e-009,   3.383409575963318e-009,   2.493865395714628e-009,       -339.2851715461339,       -15.89303469267015,         822.627341947264, /* 61 OXY */
  1.512989060582066e-009,   3.375645406087287e-009,   2.510638375996346e-009,        -197.384389531912,        612.3161967332029,        1970.683410353866, /* 61 HY2 */
  5.310257832506123e-010,  -5.663840565768641e-010,   2.172497987433713e-009,        170.2637565598232,        1432.945115000079,       -318.1623326794446, /* 62 HY1 */
  4.412828157248899e-010,  -5.319653062642984e-010,    2.20009417207216e-009,       -290.9293130666327,       -481.5619321933085,        612.8471318919909, /* 62 OXY */
  4.470055414450391e-010,  -4.335205383271641e-010,    2.21670374991566e-009,       -2257.752797358006,         -318.89334740025,        385.7688331205168, /* 62 HY2 */
  1.083834123047469e-009,   2.180712364056344e-009,   -3.85937348331243e-010,       -370.6442159261976,       -754.2431417840536,        1842.004573417689, /* 63 HY1 */
  1.178783795300589e-009,   2.203712137720018e-009,  -3.645932589750518e-010,       -198.9962469903135,       -282.0053281140746,        590.8472064372156, /* 63 OXY */
  1.195654090397999e-009,   2.299714786352672e-009,  -3.869289836967268e-010,       -1411.034818797874,         227.545124940494,        1829.012171429156, /* 63 HY2 */
 -5.419109744088746e-010,   7.434966913914375e-010,   2.140096835422582e-009,        481.2635272654404,        1312.624656302193,        1895.280503694086, /* 64 HY1 */
 -5.372925713006668e-010,   6.520806524487597e-010,   2.180368142537716e-009,        852.6994426677109,        192.4013422416984,       -641.6281202822393, /* 64 OXY */
 -6.236437463496457e-010,   6.042338945335957e-010,   2.164425673322781e-009,        1098.382786933142,        221.9327507373665,       -2095.070435352584, /* 64 HY2 */
  1.092890508221068e-009,   1.636461947195701e-009,   7.607233694577235e-010,       -954.1734182041958,        1090.106633451749,       -1086.556923019253, /* 65 HY1 */
  1.067847184011461e-009,   1.542156788610392e-009,   7.388288729635322e-010,        572.9224628180971,        451.8989852898076,        -126.165100048794, /* 65 OXY */
  1.025449251698506e-009,    1.49952077095294e-009,   8.187324732197352e-010,       -3587.775241976074,         2287.47574701608,       -1285.525630377101, /* 65 HY2 */
 -3.961311674583692e-010,  -1.025100738087201e-009,   2.842884597031612e-009,       -1249.277758943936,        -920.251105704019,        30.40196346184088, /* 66 HY1 */
 -4.477674735282132e-010,  -9.423427809611352e-010,   2.820866140505669e-009,       -549.7370801245502,       -214.6404010150235,        1019.651116425774, /* 66 OXY */
 -5.368816039267621e-010,   -9.67682962109373e-010,   2.783229097812493e-009,        -1036.75890241391,        993.0366107055164,        1347.704221409851, /* 66 HY2 */
    1.0562498110776e-009,   2.457867539072426e-009,  -4.425951073606194e-011,        388.1975288298499,       -259.8432666663832,        205.0714666549344, /* 67 HY1 */
   1.06720471745396e-009,   2.379959945137442e-009,   1.746890114653377e-011,       -792.2474523757868,       -588.2729982197776,         6.29268433398509, /* 67 OXY */
  1.042734242505153e-009,   2.407415572439656e-009,   1.104602146399038e-010,        2378.493587933792,         364.421731423883,        589.7678827468022, /* 67 HY2 */
  1.039350275643062e-009,   2.829570167504382e-009,  -3.956364020759306e-010,       -2314.637232442015,        1151.593798139567,        2057.423956051948, /* 68 HY1 */
   1.09300501683139e-009,   2.912897130933322e-009,    -3.8230276364803e-010,        168.4461275658867,       -27.64917121189056,       -317.5645086623305, /* 68 OXY */
  1.101220204953763e-009,    2.93195475134722e-009,  -2.844798737486104e-010,        1992.743184218238,        1708.331901638382,        -792.186275886874, /* 68 HY2 */
 -3.794360588198696e-010,   1.057127050229232e-009,    2.91069729530118e-010,       -1116.926723007646,        507.6605266418559,       -858.7446129092809, /* 69 HY1 */
 -3.029168422836443e-010,   1.109459359191624e-009,   3.285682486921892e-010,       -811.0395879229181,       -176.9128238742366,        -523.057127972559, /* 69 OXY */
 -3.335231268225392e-010,    1.20205611347252e-009,   3.506837736321315e-010,       -491.9174829316207,       -404.2397258400702,        894.8550939875706, /* 69 HY2 */
   3.70126019472791e-010,   3.553615706627258e-010,   4.882139686510797e-010,        -425.951443796526,       -17.92594655396899,       -225.1969085117516, /* 70 HY1 */
  4.202220547899088e-010,   4.067890467819118e-010,    4.18603610617823e-010,        244.4310495744036,       -269.8783064063461,        68.96011659710453, /* 70 OXY */
  3.595982833358039e-010,   4.262980448231854e-010,   3.415052492207832e-010,        891.8393896206459,       -259.3997868910518,       -439.6550251160614, /* 70 HY2 */
   2.38937844460808e-009,  -2.142712679594875e-011,  -7.749034385745355e-011,       -1473.026145980289,         1222.14447397654,       -1416.877119823348, /* 71 HY1 */
  2.456214367773096e-009,   5.100265794223871e-011,  -6.055264456091213e-011,        -219.947214443448,       -135.6794276818506,       -492.1069637961377, /* 71 OXY */
  2.437844272672348e-009,   1.282025639723027e-010,  -1.214025756614423e-010,        2031.122724207954,        360.3229877862112,       -564.2625781344055, /* 71 HY2 */
 -7.726192219431598e-010,  -1.205082410945604e-009,    -1.3003104418117e-009,       -897.5071072176065,         998.396793099971,        300.3220326779316, /* 72 HY1 */
 -8.623471316604582e-010,    -1.2345282801618e-009,  -1.333201823720454e-009,       -494.5734041252788,        300.8713302723131,       -181.1209971839287, /* 72 OXY */
 -8.861819610943877e-010,  -1.322068780655556e-009,  -1.291147553079312e-009,       -528.1421321796076,        658.6058351998863,        548.4480858319791, /* 72 HY2 */
  7.195057908067273e-010,  -7.437675295432578e-010,   3.242512169185732e-010,       -1098.332009976026,       -302.3479890021424,       -363.4907876514323, /* 73 HY1 */
  7.043019870032797e-010,  -6.808022360091284e-010,   2.480658104415556e-010,       -390.8175802398361,       -124.3801103656074,         -359.34534103557, /* 73 OXY */
  7.917052222257228e-010,  -6.417698036865361e-010,   2.191323747721924e-010,       -124.4813325028672,       -28.08636298095913,          567.01002704901, /* 73 HY2 */
  8.363408463138927e-010,  -1.210004944824415e-009,   9.878255153165812e-010,        234.9269290230594,       -76.66378595176664,        289.0151570564338, /* 74 HY1 */
  7.798034628598921e-010,  -1.147722584034898e-009,   1.041903531244132e-009,        387.5618655333252,       -82.26686466023469,        455.2804312037539, /* 74 OXY */
  6.834916327159181e-010,  -1.172478024644526e-009,   1.031358412252781e-009,         408.549194148813,       -532.4472700799703,        1298.746657926041, /* 74 HY2 */
  5.943641658589418e-010,  -1.129707189270516e-009,   6.607840628249484e-010,       -559.5280798822096,       -39.41811381172302,       -94.29333322729254, /* 75 HY1 */
  5.090691767129669e-010,  -1.127951078616034e-009,   7.129537944895572e-010,       -238.9854825828751,        37.30620910005604,         429.028621806659, /* 75 OXY */
  5.066743840351295e-010,   -1.04597130912754e-009,   7.701690229051087e-010,       -326.2749656643422,        628.8121165877868,       -417.4601369574743, /* 75 HY2 */
  5.147484908748284e-011,  -1.716428850168538e-009,    3.43529996580848e-009,        287.9832416864313,        175.3936478642598,        703.7078507324695, /* 76 HY1 */
  8.077931976741443e-011,  -1.639397011898553e-009,    3.37866654907493e-009,        -119.714358855169,       -611.0547749360524,       -587.7929540710101, /* 76 OXY */
  2.301175602201804e-011,  -1.560054444248849e-009,    3.39784115565952e-009,        1979.368376465716,         559.428426464237,        1001.020884820513, /* 76 HY2 */
  9.345607616830334e-010,    1.43450554115864e-009,  -8.221375887166338e-010,         455.594652422572,       -2996.423899567497,        -88.0638397926576, /* 77 HY1 */
  8.983647187074575e-010,   1.425494822404648e-009,   -9.14920418057245e-010,        130.0499844079494,        232.9747328237904,       -303.2016400905764, /* 77 OXY */
  8.239068289623473e-010,   1.490923882771143e-009,  -9.281514691970048e-010,        708.1850086518656,        1364.319553422036,        1914.534547851584, /* 77 HY2 */
  -5.86075751460003e-010,  -1.274922767964798e-009,   2.445007150809733e-009,         1687.45548000025,        478.4277870757794,        527.5613683388007, /* 78 HY1 */
 -6.029215157056636e-010,  -1.176403876034388e-009,   2.441805876097061e-009,        24.24272546461324,        174.8027731898917,       -347.3000757224618, /* 78 OXY */
 -5.610352177422254e-010,    -1.1380843301785e-009,   2.359482425924474e-009,       -168.4722417579139,         235.607588314059,        -417.189442940157, /* 78 HY2 */
 -6.540300236034923e-011,  -1.995938802977941e-010,   1.174837746299384e-009,       -422.5045915897332,       -602.9727408857832,        1214.876973503782, /* 79 HY1 */
 -1.440894818146092e-010,  -2.470482898958562e-010,   1.135385028258743e-009,         44.5030027431916,        -302.773547119959,       -90.38003830609708, /* 79 OXY */
  -2.25695119504332e-010,  -1.900795791293784e-010,   1.145136227716296e-009,        198.7531256584615,        254.5152259987984,       -1957.365127128157, /* 79 HY2 */
 -2.888226721058458e-010,  -9.647047755573652e-010,  -7.174888019920275e-011,        557.2597216343279,        -370.572331692616,        1122.351803975338, /* 80 HY1 */
 -2.133089918395679e-010,  -8.997018456834394e-010,  -8.025200906144453e-011,        378.8500771557867,        -234.698255369958,        566.0754426026778, /* 80 OXY */
  -1.56323629430307e-010,  -9.038894791897566e-010,   1.815842114895978e-012,        1390.353593045613,        -634.894409104373,       -151.5284120793542, /* 80 HY2 */
 -2.619395342465156e-009,   5.090154642201591e-010,   1.051005571453082e-009,        48.99043887232012,        65.38525227354111,       -335.8844407537263, /* 81 HY1 */
 -2.532082730176008e-009,   5.571897277875562e-010,   1.058472034135254e-009,        293.4555948363227,       -409.5399071593364,       -119.2707447359868, /* 81 OXY */
  -2.47076146472237e-009,   5.265205142525582e-010,   9.856771047492704e-010,        479.5227120898071,       -1096.932444146553,         324.655910244032, /* 81 HY2 */
  8.919173523522911e-010,   9.247386201445268e-010,    4.29702100728306e-010,        3380.204735280517,        -1157.97431002032,        659.6449578021527, /* 82 HY1 */
  9.097923023478394e-010,   9.143058013661244e-010,   3.318673295321956e-010,         359.754954956946,        245.0397318914333,       -71.52990615695141, /* 82 OXY */
  8.646983990728167e-010,   8.315695644462052e-010,   2.983821837544607e-010,        744.1196583871142,       -203.4368617233674,        513.8010084474239, /* 82 HY2 */
  3.477803461811043e-010,  -3.913018936480964e-010,  -2.444509480611263e-010,         427.653592540674,        403.6452158808247,        722.3664362539249, /* 83 HY1 */
  2.965035299705659e-010,  -3.157674740081188e-010,  -2.852582950584898e-010,       -31.94166535384229,        -312.355163706403,        -33.3761353124108, /* 83 OXY */
  1.985013803098454e-010,  -3.343762189280536e-010,  -2.782373781394899e-010,        54.83332322719036,         -377.55283314943,        1047.062622353047, /* 83 HY2 */
 -6.801719341870179e-010,  -4.400613932777652e-010,   1.485895911254245e-009,        713.0857517730511,        -953.436667385164,       -291.3600592278537, /* 84 HY1 */
 -6.075211631960592e-010,  -3.736283412516401e-010,   1.503460508400594e-009,       -9.982800030151626,       -370.1540605230406,        514.8351086400011, /* 84 OXY */
 -5.452014598172845e-010,  -3.709427317515289e-010,   1.425300147139948e-009,       -1317.681248777361,        2936.027742413164,       -457.6904114359808, /* 84 HY2 */
  1.929413216575867e-009,  -1.278428961010522e-009,    2.47364470698311e-009,       -749.5523462282041,       -207.8582804080588,       -441.7663258480522, /* 85 HY1 */
  1.988125853537445e-009,  -1.289263061172548e-009,    2.55386595438491e-009,       -478.0784315138851,       -59.72551135396915,       -620.0508565802238, /* 85 OXY */
  2.053487027115465e-009,  -1.213681640992155e-009,    2.55778589214227e-009,        314.6158333380982,       -685.8138893962845,       -1634.914858391434, /* 85 HY2 */
  3.290194875619308e-010,   -1.42767269131166e-011,    1.60458146353119e-009,        542.0464386227923,        704.5147956965316,       -150.9804497898667, /* 86 HY1 */
  3.877682879767926e-010,  -8.032378138126025e-011,   1.557822880940424e-009,        102.9563796273094,       -317.9402138692584,        730.7869851800534, /* 86 OXY */
  3.418429718920334e-010,  -1.690807163975083e-010,   1.554207382490246e-009,        149.9288283060479,       -464.1766807366064,        3274.902810131933, /* 86 HY2 */
  1.750380700063518e-009,   4.911611520459172e-010,   1.300484870915821e-009,        987.9213203459546,        328.2465540966745,        2084.453085938132, /* 87 HY1 */
  1.785627513285569e-009,    4.04204423853966e-010,   1.265898897171522e-009,       -411.7036052618667,        576.9624750782637,        -13.6683187300668, /* 87 OXY */
  1.885174290048187e-009,   4.025021442574923e-010,   1.275255257593738e-009,       -310.2407443804597,       -737.0950534758036,       -1245.311485887866, /* 87 HY2 */
 -1.008488109458454e-009,  -1.545522530117677e-009,  -1.429800833082561e-009,       -2377.123255038158,       -1178.457486099638,        1289.976308149841, /* 88 HY1 */
 -1.048123253267306e-009,  -1.575997332687238e-009,   -1.51640534672707e-009,       -273.4989423539638,        -167.995308078324,       -49.22665089241828, /* 88 OXY */
 -1.006508913931264e-009,  -1.662833787152028e-009,  -1.543380834444822e-009,        173.1771732960732,     0.004570332116623394,        96.73384868418014, /* 88 HY2 */
  2.780610671981572e-010,   1.988367790820208e-009,   5.432067099013828e-010,       -302.1852317633958,        631.6037968081463,       -383.6969023629718, /* 89 HY1 */
  2.250344200903223e-010,   1.915504015586319e-009,   4.998587577563675e-010,        -102.569796506592,        57.12388960806038,        332.6708153101996, /* 89 OXY */
  1.599124951079374e-010,    1.95551486243391e-009,   4.353741315177988e-010,       -1252.693170789576,       -659.5193144837947,        1040.442940538801, /* 89 HY2 */
 -2.259261024333882e-009,  -1.748596475990296e-009,   2.176701609983048e-009,       -2031.911307792267,       -736.4202850037522,        1456.594823459333, /* 90 HY1 */
 -2.327343345241604e-009,  -1.682370853377466e-009,   2.145413201756284e-009,       -124.3089210601069,        700.7378148619936,        291.2134254939026, /* 90 OXY */
 -2.316765411665937e-009,  -1.667627251877488e-009,   2.047073313938164e-009,        446.9789281870517,       -1186.710963479814,        59.66588193392501, /* 90 HY2 */
 -1.032167653124891e-009,   2.033500579156497e-009,  -7.264090100993455e-011,       -845.5334186619853,        -1047.14259552959,        600.4690416513185, /* 91 HY1 */
 -1.063074224735459e-009,   2.128442641901623e-009,  -7.818964437949871e-011,        439.6399322779989,       -667.7499255995139,       -174.2910970045947, /* 91 OXY */
 -1.007733444302099e-009,   2.177636736653188e-009,  -1.454007971936022e-010,        1713.821036370086,        -1726.27141282602,        89.61981782221704, /* 91 HY2 */
  6.547684413673084e-010,  -2.830177259201332e-010,   2.150650849519798e-009,        303.5540163615918,        738.9947067455611,        1210.133058160136, /* 92 HY1 */
  7.411345822828024e-010,  -2.820875961620557e-010,   2.100252193168678e-009,       -27.51258964578368,       -188.0695234695844,        619.1487693939747, /* 92 OXY */
  7.523364257432138e-010,  -1.936466769178838e-010,   2.054945212400423e-009,        86.81587506581363,       -707.6817205892291,       -373.8900185141932, /* 92 HY2 */
 -8.680543423458229e-010,   6.704902029931387e-010,   1.220806824647422e-009,        1528.929187669848,        565.3828041315487,        758.1030250360197, /* 93 HY1 */
 -8.931583731304025e-010,   6.112282618212769e-010,   1.297343155673353e-009,        616.1636506048827,         58.4156501792324,        71.27262780396764, /* 93 OXY */
 -8.984487511046969e-010,   6.654914376419744e-010,   1.381173458373168e-009,        972.6612485541794,       -742.0740490326144,        615.0977588137216, /* 93 HY2 */
  6.929103430352499e-010,  -1.018004317653676e-009,   1.833666697535374e-009,        1808.797495868942,       -251.5585106499283,          715.26591992025, /* 94 HY1 */
  7.014309776137141e-010,  -1.116594807023794e-009,   1.819268288173885e-009,        427.9590728006148,       -324.2460620081346,        360.4495334857739, /* 94 OXY */
  6.177234529961377e-010,  -1.162077294339983e-009,   1.849672131186422e-009,       -274.2166977990635,        751.6726219559713,        51.10210042018397, /* 94 HY2 */
   3.34458815231527e-009,   2.793045958512838e-009,   7.671888862179262e-010,        934.2476546005376,        1579.882918832721,           336.3222458249, /* 95 HY1 */
  3.397445176700214e-009,   2.771167022681828e-009,   6.851678697643812e-010,        300.1866765718372,        965.7297715737132,        88.97516756945303, /* 95 OXY */
  3.337020490367734e-009,   2.771547986141626e-009,   6.054890636208087e-010,       -779.9207196563628,       -1952.195366672141,        861.8777294870684, /* 95 HY2 */
 -2.979568408317928e-009,   2.974878862119825e-010,   1.412650030905928e-009,       -30.36566622069668,       -1100.709175219643,       -700.2781980600486, /* 96 HY1 */
 -2.986892224054564e-009,   2.233356465972615e-010,   1.479341912594619e-009,        15.85501683396096,       -76.61722337140672,        452.3672918540437, /* 96 OXY */
 -3.054615412706677e-009,   1.568085719340685e-010,    1.44791393054757e-009,        2969.863992493352,       -2679.897767256207,       -533.5401393078118, /* 96 HY2 */
  1.850008014354848e-009,   1.541366082020075e-009,  -9.808021592022018e-010,        -278.821304958107,        439.8574525780318,        23.01572935829064, /* 97 HY1 */
  1.752698225461219e-009,   1.526434737151358e-009,  -9.983480967916104e-010,       -47.96505935408371,       -524.5593651979819,       -453.8601070232782, /* 97 OXY */
  1.701425021692108e-009,   1.537278485763988e-009,  -9.131807809287844e-010,       -698.3244499518703,         61.7115567159297,        -917.160512871002, /* 97 HY2 */
   1.04911021572935e-009,   2.952389943946608e-009,   9.389722196534254e-010,       -454.3131069957452,        -77.1668142299565,        1331.434027623875, /* 98 HY1 */
  9.631835138147048e-010,   2.963095811559349e-009,   8.889523606298989e-010,        493.7075646893375,        372.9653583424415,       -218.2901728570266, /* 98 OXY */
  9.717971247622014e-010,   3.038309941130366e-009,   8.236175275060568e-010,        1512.758235193846,          1441.2049402563,        1130.526158835719, /* 98 HY2 */
   3.57392151215682e-010,    -9.5208025525675e-010,   1.270915903636797e-009,       -2016.914913015204,        2496.380928861732,       -1061.336750534047, /* 99 HY1 */
  3.243118388647794e-010,  -9.452778411516454e-010,   1.365040396183156e-009,        326.9728826898832,       -214.7083803295246,       -4.559302736008254, /* 99 OXY */
  2.710147862126152e-010,   -1.02682439473599e-009,   1.387614358553731e-009,        60.63684579484179,       -423.5685393701252,       -1366.069719537832, /* 99 HY2 */
  3.750167748455948e-009,    3.72069872137146e-009,  -6.320252096787605e-010,        295.6464338812231,       -1281.346319511807,       -688.0981477584976, /* 100 HY1 */
  3.738838195518419e-009,   3.756735683476542e-009,  -7.246155903479972e-010,        474.0705168400596,        312.5471079842157,       -97.46102485081228, /* 100 OXY */
  3.818185693962016e-009,   3.731675618654746e-009,  -7.800769971267038e-010,        499.0146130705568,        1113.330125447271,        -426.988162664355, /* 100 HY2 */
  2.275218164493938e-009,  -1.145935044384396e-009,   5.477632745492209e-010,       -142.4217049850794,       -325.6479210779824,        382.8150718166278, /* 101 HY1 */
  2.259822442481947e-009,  -1.149327993323244e-009,   4.490137953544954e-010,        -185.845508384145,        408.4249098583195,         362.993031152037, /* 101 OXY */
  2.329155321016476e-009,  -1.207475861425504e-009,   4.064488169182312e-010,       -749.2624768323203,       -9.362066955942644,        12.38493897046784, /* 101 HY2 */
 -1.235737252255844e-010,   4.408030007448316e-010,   9.475559488417484e-010,       -24.16738212461282,       -2051.550337151978,       -419.4617539710983, /* 102 HY1 */
 -1.858838188172822e-010,   4.345285363579498e-010,   1.025518012265268e-009,        225.4351702153558,        33.93449254785736,        -37.5142340858903, /* 102 OXY */
 -1.340268011606033e-010,   4.440819219896676e-010,   1.110486138893598e-009,        730.2348373685563,         657.976817761399,       -413.4519924198759, /* 102 HY2 */
  3.856486182300965e-010,   1.302759041354617e-009,   1.749413442222658e-009,         174.360309711822,       -446.7730811197938,       -318.5047368020003, /* 103 HY1 */
  3.909057960441246e-010,   1.280614367014232e-009,   1.652038006512247e-009,        375.5964918069736,        -7.03552700380846,       -408.2644483379742, /* 103 OXY */
     3.814767786051e-010,   1.181819066213919e-009,   1.639766813827273e-009,       -1119.575835018936,        201.2821614700629,       -989.8791274620724, /* 103 HY2 */
  1.419498294860993e-009,   1.970712245866952e-009,   2.620253060509186e-009,        -1808.58343443436,        3047.740308414051,        1313.143545660268, /* 104 HY1 */
  1.423645527314278e-009,   1.948772550923964e-009,   2.717728443292678e-009,       -668.8617230702043,        522.9254917564091,        716.9613720735123, /* 104 OXY */
   1.51894896063319e-009,   1.935236840715973e-009,   2.744821618573432e-009,       -338.2422534559222,       -57.95218175945077,       -713.2405377361052, /* 104 HY2 */
   -3.4394391594624e-010,   2.267960748123408e-009,   1.868851363104146e-009,        798.9797889745861,       -230.2566021980854,        -3643.34358750071, /* 105 HY1 */
 -2.622508938258912e-010,   2.214274935076554e-009,   1.889924654597015e-009,        270.4076018368374,        142.3832980437711,       -524.5925521113963, /* 105 OXY */
 -2.864234063023226e-010,   2.117369469397712e-009,   1.894926690396985e-009,        62.01578301803833,        271.3943222633688,        1103.125460222292, /* 105 HY2 */
  8.037875738173886e-010,   5.043378271388718e-010,  -5.126065728294908e-010,        47.01371521816591,       -589.7676500548472,       -803.4295751147645, /* 106 HY1 */
  7.141199975743249e-010,   4.667260348476062e-010,  -4.892596490207075e-010,        113.0937896227438,       -505.2173365922038,       -411.6615463638558, /* 106 OXY */
  6.598942398705229e-010,   4.556802375832259e-010,  -5.725516508230582e-010,        -1054.60692459498,        1576.681874927596,        54.70241879032802, /* 106 HY2 */
  2.234371810506896e-009,   2.686992030169198e-009,   9.668247692797398e-010,        110.1539817475733,        783.1178378461556,        744.2701908032041, /* 107 HY1 */
  2.200637608181253e-009,   2.694364840921182e-009,   1.060673824975938e-009,       -1.987169156305882,        562.9827402635482,        721.4187244640152, /* 107 OXY */
  2.118883278332646e-009,   2.751929130785766e-009,   1.062280711216927e-009,        -430.584182928556,       -41.42300597563812,        653.5409944101083, /* 107 HY2 */
 -1.386003419914595e-009,   9.838060750143928e-011,  -1.058903972096007e-009,       -133.3515839217412,        1101.397292193047,       -1259.478096750276, /* 108 HY1 */
 -1.441873110592295e-009,   1.682249753424657e-010,  -1.014178789654058e-009,       -268.4122086497732,        56.49500594297962,        222.0358315260949, /* 108 OXY */
 -1.408684680282727e-009,   1.828162905547152e-010,  -9.209821072347582e-010,        1510.467190694206,        188.2035372129703,       -422.4157186511962, /* 108 HY2 */
  1.280880300790909e-010,   1.395162542068222e-009,   3.224181884987702e-009,        54.68881962579844,       -780.7618348782415,        -910.496616435803, /* 109 HY1 */
   1.08666512864449e-010,   1.428141852563967e-009,   3.131795902564462e-009,        69.22252422238783,        480.4955217778945,        -468.151709524586, /* 109 OXY */
  1.944333837984684e-010,   1.438801003029241e-009,   3.081492559121438e-009,       -182.5939449564336,       -2903.880755471887,       -1679.166227291846, /* 109 HY2 */
  1.635839912268506e-009,   1.326176550712234e-009,   2.974700227024088e-009,        272.2111888501994,        1225.374476165728,       -229.9201476047592, /* 110 HY1 */
  1.674307135549953e-009,   1.365505023494921e-009,   3.058207974937669e-009,       -197.7202890998685,        236.8097444856124,        457.1206586638714, /* 110 OXY */
  1.686056315787512e-009,   1.464093351841941e-009,   3.046279090769826e-009,        1127.344752948072,          175.47844792795,        1206.689909042577, /* 110 HY2 */
  1.213336869673889e-009,  -1.179108185643996e-009,   1.406769829743786e-009,        2050.722231861826,        1281.877652839773,        920.4736737517308, /* 111 HY1 */
   1.24037732140111e-009,  -1.249738430563413e-009,     1.4721929332548e-009,       -73.64639733385889,       -632.6413505382004,       -232.0609722219594, /* 111 OXY */
  1.216567282816545e-009,  -1.220310535174296e-009,   1.564751460876889e-009,        -86.2981263662288,       -1875.821793627823,         164.538715172291, /* 111 HY2 */
  6.954602458676299e-011,  -1.637951387464164e-009,   1.206044187102008e-009,        610.8936503091091,         1728.86258885988,          35.802586010031, /* 112 HY1 */
  1.452651666785382e-010,  -1.584236958852742e-009,   1.243210995948541e-009,        1039.793806183001,        862.9008647160309,         420.802139022198, /* 112 OXY */
  1.545984333840562e-010,  -1.602998059131111e-009,    1.34099090824626e-009,        1222.579937183012,       -299.4927618680654,        184.0089738520199, /* 112 HY2 */
  -1.38842540101053e-009,   -1.11062117399682e-009,   2.138386114898251e-009,       -549.8095163261386,       -290.6740939162978,       -862.4849165884938, /* 113 HY1 */
 -1.368208590733676e-009,  -1.094142052382311e-009,   2.041847422484049e-009,        42.91081687215231,        805.0520752896513,        -555.582182034217, /* 113 OXY */
 -1.373259668740288e-009,    -1.1804164509279e-009,   1.991536243144184e-009,        1755.976802074342,        1137.788920041819,        -1316.15761567588, /* 113 HY2 */
 -8.921251406313811e-010,   3.699957112260669e-010,   2.510378452143974e-009,        161.9367041725568,        -2031.67953551966,        494.1804646917399, /* 114 HY1 */
 -9.654498592532958e-010,   3.760686345608474e-010,   2.442653971775438e-009,        572.7283858361241,        142.5893487683758,        226.0493776054124, /* 114 OXY */
 -1.053800731032756e-009,   3.645784090175974e-010,   2.488064302899288e-009,        246.5750923840992,       -432.2743521762037,       -548.2708378594209, /* 114 HY2 */
 -5.284365255835126e-010,  -5.657353594024184e-011,   1.493665577847593e-009,        489.2733337389978,       -1775.273073468344,       -445.7932753883077, /* 115 HY1 */
 -4.913445759360522e-010,   -1.98934974565277e-011,    1.40834991663475e-009,       -33.05151273734652,       -388.7230962254013,       -83.57346297575708, /* 115 OXY */
 -5.630980014202623e-010,   2.964066025431986e-011,    1.35938244629821e-009,       -1085.714882639079,       -1995.868542976759,       -185.7179244737036, /* 115 HY2 */
 -1.155342638149152e-009,  -1.194869570009669e-009,   5.290748097584416e-010,       -922.0880002838766,        56.45214042393023,        905.5098580022528, /* 116 HY1 */
 -1.111992592625018e-009,  -1.112800553351192e-009,   5.662964440289487e-010,       -278.4792615839057,       -73.08784580402619,        445.8664823804689, /* 116 OXY */
 -1.012761329461931e-009,   -1.12512283144009e-009,   5.674444308578361e-010,        -234.168488938173,        18.54413535989029,       -1539.833822623134, /* 116 HY2 */
  2.023931704903489e-009,   8.437113038972859e-010,   5.622359079687613e-010,            1231.32452688,       -2235.083034895049,        134.6226227736038, /* 117 HY1 */
  2.020616712927894e-009,   8.468785883548204e-010,   6.621307483628832e-010,        551.2496195593171,         105.285005739495,        52.73229159694966, /* 117 OXY */
  1.925537801014391e-009,   8.408964955462038e-010,   6.925316541877949e-010,         201.065559591973,        2628.675418466134,       -490.1381985059554, /* 117 HY2 */
  1.799388507999348e-009,  -1.396803940577669e-009,   2.562903149175169e-010,        1979.922632691809,        469.1623730240962,        645.8081378273694, /* 118 HY1 */
  1.837846953479326e-009,  -1.474646669065999e-009,    2.06677236069728e-010,        599.8028081079007,        330.4654645438474,        -219.872055339059, /* 118 OXY */
  1.778036977766248e-009,  -1.498189917073462e-009,   1.300712162118762e-010,        403.0416718279872,        2528.022908039739,       -758.4573144764269, /* 118 HY2 */
  1.535085796721065e-009,   1.415934259218598e-009,  -1.353896728774753e-009,       -1514.060362268112,        -35.9933980559459,        493.8457719459691, /* 119 HY1 */
  1.522375536757689e-009,   1.324460602152391e-009,  -1.392250598913216e-009,        -318.421820822067,       -106.2969768384662,        255.5707695001414, /* 119 OXY */
    1.5642631786674e-009,    1.25680693398133e-009,  -1.331683226192876e-009,        1480.481068200308,        366.4071182082144,       -444.2174345880778, /* 119 HY2 */
 -3.757557587723064e-010,   3.276575652708314e-009,  -5.971837030031063e-010,         480.295797383037,       -209.7573683642842,        1049.815758445227, /* 120 HY1 */
 -4.683229234421658e-010,   3.274229824759067e-009,  -5.594237094091473e-010,        66.34036165144163,        -454.456860313509,        28.25838727961816, /* 120 OXY */
 -5.006644255696681e-010,   3.179655295417049e-009,  -5.563115333420278e-010,        1209.572415708586,       -796.9324771937739,        1895.914079353008, /* 120 HY2 */
 -1.831416643823976e-010,  -2.803405913589508e-009,   1.290076148110913e-009,          391.02798526402,       -1120.196589280135,         1872.91806298981, /* 121 HY1 */
 -9.659195324267174e-011,  -2.758215150765694e-009,   1.311684999469375e-009,        143.4953212373683,        532.5651006818888,       -494.8994909979007, /* 121 OXY */
 -1.149704211030094e-010,  -2.669381629713209e-009,   1.353766319036674e-009,       -448.0177712441988,       -659.5006840305158,        1805.161312381976, /* 121 HY2 */
 -4.491263660825751e-010,  -1.086439813139381e-009,   1.291362998948622e-009,       -186.8807723823108,        1826.282877005883,        798.2864101323737, /* 122 HY1 */
 -3.513557295647927e-010,  -1.080290454211091e-009,   1.311440051015094e-009,         228.163324546499,       -691.2146352981971,       -354.2014530457285, /* 122 OXY */
 -3.322861259578024e-010,  -1.125663468957954e-009,    1.39848969059537e-009,        274.2862290915433,       -2136.561543144035,       -1110.019494437855, /* 122 HY2 */
  2.351850951515664e-009,    5.77596832712928e-010,  -1.044910832780964e-009,       -944.8484340322196,       -1929.123208163986,       -959.4141318842292, /* 123 HY1 */
  2.378532035510798e-009,   5.156535720371742e-010,  -9.710789141163308e-010,        323.2119460303664,        73.18967372680083,        286.5097862980646, /* 123 OXY */
  2.297058729524532e-009,   4.864741117933223e-010,  -9.209724323844448e-010,        1151.324810763564,        1209.327890132282,        2325.253139690642, /* 123 HY2 */
  8.869955766297609e-010,   4.898164323071352e-011,    5.17959299596511e-010,        301.9760335562795,        40.04042668470896,       -45.09636059690479, /* 124 HY1 */
  8.627868081491926e-010,   5.605674322403775e-011,   6.147264349717977e-010,         205.333535877023,       -129.0390724181616,       -56.81340946555697, /* 124 OXY */
  9.263966587738808e-010,   1.181403998008437e-010,   6.605457157339207e-010,       -1030.142180883345,        1327.896033076791,       -295.5103874596136, /* 124 HY2 */
 -7.938371450724889e-010,    1.92319575855177e-009,   1.595823449895891e-009,       -958.4353439726374,       -30.52112316846864,       -125.9827568982273, /* 125 HY1 */
 -8.442224290231895e-010,   1.937076481894569e-009,   1.510567079572683e-009,       -540.8331947686594,       -470.8936381815955,       -445.8574083319855, /* 125 OXY */
 -9.392771185580738e-010,    1.96011277018482e-009,   1.531398196183062e-009,       -828.8802885312072,       -1109.068058097562,       -1044.334875306043, /* 125 HY2 */
  5.753645075024436e-010,   4.213051637587256e-010,   2.617057281619878e-009,       -1040.006242928154,       -2280.457707386006,        1444.560402155768, /* 126 HY1 */
  6.260947279026832e-010,   3.901492690322998e-010,   2.536709568759214e-009,       -31.44416772082646,        239.7889014944131,        1080.750903756435, /* 126 OXY */
  5.617887212793438e-010,   3.599431317338616e-010,   2.466336798587314e-009,        911.2528995316438,       -525.9176322607605,        541.7152370397392, /* 126 HY2 */
  1.067669460165724e-009,  -1.797200814936764e-011,   3.071330732111713e-009,        1428.934456407896,       -1484.851987164839,        734.0990266457377, /* 127 HY1 */
  1.121112634922321e-009,   4.885706646833808e-011,   3.123077247527516e-009,        147.0742996139152,        83.25325165466852,        54.87494471178322, /* 127 OXY */
  1.077902631826295e-009,   6.451934619454021e-011,   3.211889344824127e-009,        643.6642063477057,        -959.504162674128,        484.6505472054275, /* 127 HY2 */
  5.465027555514914e-010,   1.410780849315595e-009,   2.458021469589836e-009,       -149.2882680705184,        315.1845765536408,       -217.1511168317041, /* 128 HY1 */
  6.266230023135683e-010,   1.351286684167298e-009,    2.46443944217411e-009,       -393.8078194729615,       -38.82521425215494,       -437.1845127002731, /* 128 OXY */
  6.009440411351741e-010,   1.257968088398316e-009,    2.43929524529878e-009,       -860.5747439519409,        86.72914036811859,       -428.7886497810449, /* 128 HY2 */
  1.753617924010121e-009,   2.500305657042368e-009,   7.425656181154734e-010,         1155.40986443797,       -1849.191758628023,       -640.0718672299637, /* 129 HY1 */
  1.670237606853703e-009,    2.45123896344356e-009,   7.678674419662543e-010,        154.8183758182636,       -13.10040640214675,        -332.662209754688, /* 129 OXY */
  1.607665663382188e-009,   2.513719860426137e-009,   8.145674355685031e-010,        2212.093767772266,        1080.616099277653,         999.069406167131, /* 129 HY2 */
  1.466355765473446e-009,   2.474878783377564e-009,   1.559830042746798e-009,        1611.776887918261,       -2309.874668235772,        409.5208726957435, /* 130 HY1 */
  1.540579135850731e-009,   2.541856703176848e-009,   1.557627873544131e-009,        -184.187945491958,        -303.218914412674,        72.40514500004873, /* 130 OXY */
  1.561485472283018e-009,   2.571517220114401e-009,   1.650811449521754e-009,       -191.1420548678996,       -153.5990927385376,        26.40692387575646, /* 130 HY2 */
 -1.252363644607455e-009,   1.101375345168102e-009,   1.777649369780988e-009,       -731.5258967021721,        317.2047832560491,       -673.7591228361282, /* 131 HY1 */
 -1.184667078255884e-009,   1.077270211890154e-009,   1.708107170693958e-009,       0.9869249036142946,       -35.90058911582766,        156.8511066919274, /* 131 OXY */
 -1.120120725773648e-009,   1.152811647138175e-009,    1.69682621586255e-009,        1011.770653832099,       -686.2185508279974,        1512.737601898076, /* 131 HY2 */
 -1.217682708603058e-009,  -1.180850735277528e-010,   1.401481829710825e-009,       -2801.823041758046,        639.9667304796351,       -735.6104080240281, /* 132 HY1 */
 -1.214602530008133e-009,  -1.331327209295898e-010,   1.500295193217462e-009,       -621.3801747260735,        230.6885381246143,       -853.4173261361743, /* 132 OXY */
 -1.306480598696816e-009,  -1.229626287525526e-010,   1.538439520700748e-009,       -82.57792745986939,       -1591.745524022185,         975.895891360052, /* 132 HY2 */
  7.957861094605614e-010,   1.076117651081847e-009,   1.428769486316675e-009,        1302.633775155242,         85.0783742659965,        1429.124026867705, /* 133 HY1 */
  7.713428944589713e-010,   9.832579025520534e-010,   1.456690738856091e-009,         200.918824633242,         95.6479485555697,        518.0974489540306, /* 133 OXY */
  7.411539669115863e-010,   9.839746801237174e-010,   1.552022343707012e-009,        1672.586810727614,       -713.2661853453274,        998.2168766779118, /* 133 HY2 */
  4.797936021693991e-010,   9.692908078809996e-010,   1.997290882895993e-009,         1976.83135729582,       -406.2179389199781,        763.6309040535514, /* 134 HY1 */
   5.10455917989096e-010,   1.055793804054404e-009,   1.957578743170243e-009,        220.3122344000824,         64.2999538977332,        410.7019749940261, /* 134 OXY */
  5.817332724747941e-010,   1.095578966222504e-009,   2.015342742415265e-009,       -1356.128370094174,        2309.401174327616,        842.9942374304366, /* 134 HY2 */
  5.723117842375845e-010,    2.32111002510421e-009,   1.911381083264556e-009,        725.9473591179554,        1075.372230544647,        650.8530658513189, /* 135 HY1 */
  5.267977987098848e-010,   2.401453931222864e-009,   1.872997051806366e-009,       -701.3532417495898,       -4.452939594520392,        59.89525146388376, /* 135 OXY */
  4.508555374288923e-010,   2.372164965488213e-009,   1.814902908919801e-009,        546.4440867430216,       -1817.491202383624,       -680.3957178718414, /* 135 HY2 */
  8.808981402064003e-010,   1.939871449902969e-009,   1.137866542258328e-009,       -2424.489791923567,       -285.1831548580055,       -1807.948606695926, /* 136 HY1 */
  8.590819908658173e-010,   2.037237758846264e-009,   1.131244041921671e-009,        -68.2627457219058,        333.7700092497198,       -737.2002799262699, /* 136 OXY */
  9.324053654239084e-010,   2.090752985714902e-009,   1.173194050701337e-009,        912.5670258851586,       -1928.284492706181,        479.1575618681494, /* 136 HY2 */
  1.007804543943707e-009,    3.23158565348515e-009,     2.1120048221877e-009,        385.6648285220047,       -2483.385378607452,        1255.123152996011, /* 137 HY1 */
  9.248414815198786e-010,   3.215081993124338e-009,    2.05866853559072e-009,        1045.577720410244,       -394.5862633681378,        -453.863961619999, /* 137 OXY */
  8.657083316013719e-010,    3.29559501400192e-009,   2.063242816942208e-009,        2355.616290698081,        581.1997629424255,       -547.3212886560083, /* 137 HY2 */
  2.053584390639228e-009,    2.71860771032822e-009,   -1.64387873048255e-009,        964.4716937997726,        531.2932149366657,         965.664679801974, /* 138 HY1 */
  2.023478616320919e-009,   2.637024750467484e-009,  -1.594504010097953e-009,       -298.2708959553406,        415.8973699472104,        17.73873725009529, /* 138 OXY */
  2.095710161359609e-009,    2.60730664679712e-009,  -1.532058528151198e-009,       -2244.506718351807,        842.2782306920943,        2512.706685632146, /* 138 HY2 */
   2.83113764708523e-009,   6.332224494621235e-010,   2.086902310187646e-009,       -1693.673657640708,       -284.2020998027598,        2075.386496757752, /* 139 HY1 */
  2.917186274521032e-009,   6.110623173259746e-010,   2.132777818469548e-009,       -307.8745591980538,        360.6832100366148,       -172.1830504202105, /* 139 OXY */
  2.903131686809556e-009,   6.110399795698682e-010,   2.231785232663058e-009,        1863.422613103845,       -1473.890606754717,        156.3047779558842, /* 139 HY2 */
  9.610958517199784e-010,   1.329597379137706e-009,   6.581847361303509e-010,       -2702.238437286304,        576.5967997359867,        3012.379898324228, /* 140 HY1 */
  8.891424638762169e-010,   1.269832528788051e-009,   6.228170582239539e-010,       -745.1633513280112,         267.080904737749,       -564.3236590454936, /* 140 OXY */
  8.048731280604825e-010,   1.322400442303505e-009,   6.111854843918543e-010,       -1888.198405873163,       -1297.542884410008,        538.7891576229082, /* 140 HY2 */
  1.803870888928654e-009,   9.400079719548712e-011,   9.087338319508702e-010,        1071.381128612076,       -885.5011971179318,         -1010.2375342048, /* 141 HY1 */
  1.783627399231598e-009,   1.253627609417374e-010,   8.159619161850283e-010,        74.34713988791405,          418.94501976061,       -360.1069411339272, /* 141 OXY */
  1.821473583834845e-009,   2.169760602704502e-010,   8.027456498028213e-010,        751.9212871260111,        224.7241642506048,        218.1635298969988, /* 141 HY2 */
  2.326721783577674e-009,   2.040474082181793e-010,  -1.691003352458482e-009,       -295.3231651755183,       -1176.423048817356,        -954.620478315934, /* 142 HY1 */
  2.243821727767708e-009,    1.54678015202283e-010,  -1.664730861523136e-009,       -482.5363596998494,       -452.3332401982104,        -173.566247402084, /* 142 OXY */
  2.164665611517044e-009,   1.962421742625356e-010,  -1.709527400917664e-009,       -667.9807419676256,        838.3236571348339,        1329.546347547915, /* 142 HY2 */
  2.277953477028558e-009,   2.092181463639309e-009,    1.33672261393324e-009,        142.6609944098782,        1485.752641387513,         91.9859456661689, /* 143 HY1 */
  2.364634561351562e-009,   2.050878628112132e-009,   1.308785803127364e-009,       -345.6086751635327,        211.3458605817903,         443.363659364029, /* 143 OXY */
  2.390202217717766e-009,   2.085015472272878e-009,   1.218337083851273e-009,       -631.7195565245237,       -1932.731048840381,       -461.9208043157774, /* 143 HY2 */
  6.847370924803547e-010,    -3.5764871150009e-010,   1.186728626389534e-009,       -838.0471770893757,        542.6183338596898,         1282.78894014014, /* 144 HY1 */
  7.329167227906047e-010,  -4.322105589118593e-010,   1.232763983819576e-009,       -498.2831679252857,        131.2958571652457,        268.1312571683099, /* 144 OXY */
  6.912130695401549e-010,  -4.486734386053408e-010,    1.32214954695907e-009,       -541.7320485778587,       -1401.078358965656,         -27.552929063603, /* 144 HY2 */
  8.275047486418211e-010,  -4.877753357810351e-010,  -2.070536884087094e-010,        688.4111440245717,        400.3998691049153,        437.1139287643706, /* 145 HY1 */
  8.675935979642398e-010,  -5.724785765380207e-010,  -2.419570573760364e-010,       -29.87555344918778,        352.5938634914584,       -279.2586460614002, /* 145 OXY */
  9.204161297672632e-010,  -5.527936216768896e-010,   -3.24554168015784e-010,        2364.591264803778,         1054.14223925231,        1391.951167117861, /* 145 HY2 */
  4.645170626103939e-010,   3.435701723034064e-010,   1.607341388622519e-009,        -942.107035426963,       -1102.696777209165,       -504.7724715301593, /* 146 HY1 */
  5.170017637163304e-010,   2.744610078142982e-010,   1.657033234366443e-009,       -354.2039404551718,        -379.544234619801,       -114.8553077961077, /* 146 OXY */
  6.086650093480053e-010,   2.675368254520228e-010,   1.617664400260247e-009,       -605.4455004586807,       -490.3265332772506,       -682.8693023138986, /* 146 HY2 */
  1.910885353339192e-009,   1.283521323694545e-009,   1.938041930143516e-009,       -67.31253080767212,       -1581.404750658255,       -855.5511561682952, /* 147 HY1 */
  1.941379751076106e-009,   1.336698492314521e-009,   1.859033826481767e-009,        166.8746890626495,       -76.05433300184548,        236.9073379020866, /* 147 OXY */
  2.034916945070307e-009,   1.368556687278146e-009,   1.874391194529419e-009,        -751.173371185405,        1746.591107049964,         2176.50864900194, /* 147 HY2 */
  2.727471049020729e-009,   3.051204847300462e-009,   2.901857865149772e-009,       -568.7385397692765,       -935.7252437437402,       -342.7903549169802, /* 148 HY1 */
   2.69224664487472e-009,    3.04889136721203e-009,   2.995420085957316e-009,         196.231931316201,       -84.73248309201752,       -29.99116581988274, /* 148 OXY */
  2.661173048672712e-009,   3.140262879905432e-009,   3.021605545418543e-009,       -1300.819428160032,       -227.9311241518838,       -1270.539786353996, /* 148 HY2 */
  1.223175219174834e-009,   1.109531627630463e-009,   1.118813804364842e-009,       -245.4279710959278,       -2698.671577165908,       -779.4187196418429, /* 149 HY1 */
  1.189805880967619e-009,   1.137685175337676e-009,   1.028847894048139e-009,        55.72515527005192,       -1141.475821539978,       -411.1846040201694, /* 149 OXY */
  1.230968382343788e-009,    1.22535853202368e-009,   1.003967560290986e-009,         995.970136395811,       -1244.620703815922,         758.182645843962, /* 149 HY2 */
 -5.277936337899503e-010,   4.844376727743144e-010,  -7.007402990192409e-010,       -3228.071117303958,       -741.2622493594444,         486.340915789334, /* 150 HY1 */
 -5.000975720085601e-010,   4.577922838821326e-010,  -6.084204668455077e-010,       -363.4544479069565,        92.63163500321096,       -107.3096029416462, /* 150 OXY */
  -5.32318017069217e-010,   5.261433689241476e-010,  -5.429225883022087e-010,       -112.7036059083283,       -359.6515859439574,        490.4124981249806, /* 150 HY2 */
 -1.898492282598283e-010,  -8.524257254962719e-010,   1.085165047448992e-009,        1302.105260469977,        1201.556549733998,       -678.4393375136538, /* 151 HY1 */
 -2.250984859058344e-010,  -8.105410522457622e-010,   1.001480138522088e-009,       -651.8036070748625,       -199.0390575738169,       -573.7281088943058, /* 151 OXY */
 -3.119496984028286e-010,  -8.534904640077014e-010,   9.767371613950498e-010,       -1251.981937636465,       -615.4679437047124,        2174.141995243224, /* 151 HY2 */
  1.094016085972238e-009,   2.205596788768094e-009,  -3.861398762831596e-011,       -2578.080894427846,        327.8788938955368,       -719.0073115121608, /* 152 HY1 */
  1.108763172393384e-009,   2.122838432303774e-009,  -9.277650120933576e-011,       -494.9713871796927,        437.7899194896106,       -340.5154671542727, /* 152 OXY */
  1.127300512108582e-009,   2.148029375386112e-009,  -1.877595703475082e-010,       -3596.220749883875,        409.6144450962702,       -979.6339307915702, /* 152 HY2 */
   2.96889728755213e-010,  -1.243347166437238e-011,   9.131713244300908e-010,        1646.615258619472,       -1348.715175445501,        545.0219524912914, /* 153 HY1 */
  2.244731838638198e-010,   -3.32533748463976e-011,   9.789161010856126e-010,        441.1217871588435,        205.2788826258244,       -273.4334998315888, /* 153 OXY */
  2.570360786324588e-010,  -1.021570909450031e-010,   1.043661261461428e-009,       -1408.773283127826,        282.6857718192141,        756.6645310380552, /* 153 HY2 */
    9.0838741474463e-010,    4.02113408365382e-010,   1.053366310488132e-009,        803.1547918930727,        776.3354355540002,        1244.307063643921, /* 154 HY1 */
  9.330217498918566e-010,   4.598200903073647e-010,     9.7550047464788e-010,        117.7114587815343,       -65.87070464175349,        397.2019804702078, /* 154 OXY */
  9.812539936725226e-010,   4.049019477768828e-010,   9.072535092876198e-010,        639.8850087667852,       -528.8208216776005,        1134.994944356382, /* 154 HY2 */
  4.710976798152529e-010,  -8.547933477207666e-011,   2.517072490599704e-009,        -1055.28899645859,        1019.885058311339,        267.4032182438203, /* 155 HY1 */
  3.955773462773441e-010,  -1.291127001192374e-010,   2.565988833991332e-009,       -455.9195915909862,        8.447832581142611,         297.617206599734, /* 155 OXY */
  3.088420501477126e-010,  -9.028181595767784e-011,   2.534857646404957e-009,        -817.245309076381,        1012.712739394254,        2508.562154285634, /* 155 HY2 */
 -5.574168658008156e-010,   1.002852122653281e-009,    4.04237866941116e-009,        577.8675594985582,        514.9918262362555,       -482.0114460718136, /* 156 HY1 */
 -5.241394384826026e-010,   9.671267081402664e-010,    3.95510721022229e-009,        378.6944867466622,        693.9096329958545,       -631.4690712002078, /* 156 OXY */
 -6.007121506166048e-010,   9.580421493228684e-010,   3.891435469649752e-009,        159.6850892208076,        392.0971250263335,       -325.9351696712094, /* 156 HY2 */
  1.418251397334716e-009,  -4.722228166762291e-010,    2.81405307812589e-009,        611.9517784902868,        499.8872609850387,        1952.680010016252, /* 157 HY1 */
  1.323376151520961e-009,  -4.587430488028721e-010,   2.842636002921306e-009,       -35.64404389819068,        140.8736504896045,        10.22511547229627, /* 157 OXY */
  1.303267188286699e-009,  -3.608557118769824e-010,   2.846337200997677e-009,       -484.9571677825429,        65.06130653539414,       -400.4734796840808, /* 157 HY2 */
  1.150569316697259e-009,   1.967538000761432e-009,   1.079639091170739e-009,        653.1239750377694,       -253.0459272360619,        708.6858011040799, /* 158 HY1 */
  1.164886138028924e-009,    1.96095575891654e-009,   9.808883812023884e-010,        -419.498028748085,       -154.2941002496418,        543.5889860132057, /* 158 OXY */
  1.262232614925838e-009,   1.947861290080602e-009,   9.621214336172324e-010,       -476.7024569160142,        983.9218233780254,       -581.4810332709957, /* 158 HY2 */
 -4.718365915613631e-010,     2.5895804416702e-009,   1.813506937238609e-009,       -941.9642020359904,        837.6382701916328,         43.2373139897553, /* 159 HY1 */
 -4.074659069818302e-010,   2.662776888218632e-009,   1.791175453198246e-009,        109.0831719955427,       -255.5464337236838,       -539.8194989300755, /* 159 OXY */
 -4.486195355395808e-010,    2.72387208947934e-009,   1.723545909066673e-009,        1074.616016505566,       -331.6187239195504,       -1201.167850076107, /* 159 HY2 */
  1.834323667139121e-009,   2.871976793252388e-009,    2.23943285388458e-009,        260.0421630015418,        1947.558677292154,       -1341.022689792324, /* 160 HY1 */
  1.857914087522038e-009,   2.904079132796804e-009,   2.331154889742563e-009,         351.494994794185,        78.18550769238298,       -699.6018357624341, /* 160 OXY */
  1.935501461383591e-009,   2.966918125598786e-009,   2.325545932539799e-009,       -338.1799150036556,        1036.148971363935,        378.7366288040353, /* 160 HY2 */
  3.656029942615474e-009,   6.714934443585173e-010,   1.671939448643405e-009,       -76.77678396395531,        1060.959089904799,          1068.3495439602, /* 161 HY1 */
  3.651449654202814e-009,   6.210723887523934e-010,    1.75817597286384e-009,        135.1344174567112,       -307.0170498118862,          287.09668248639, /* 161 OXY */
  3.644787006196796e-009,   5.231157450041749e-010,   1.739199538342742e-009,       -600.8127515239402,        34.63909825860854,       -1258.281498585996, /* 161 HY2 */
   2.79193880619222e-010,   2.785908270623956e-009,   1.739120623515096e-009,       -273.7223069647148,        1758.017007311903,        1244.313571300034, /* 162 HY1 */
  1.975222937145188e-010,   2.733319617973152e-009,   1.762873212949558e-009,        491.1841500169861,       -183.0881217522154,       -350.6599548656626, /* 162 OXY */
  1.384486722600132e-010,   2.725413243358596e-009,   1.682575040562478e-009,        1306.148036346205,        123.3049629740873,       -983.9883613257596, /* 162 HY2 */
   1.26110602437688e-009,   2.321131451124594e-010,   3.850237127602014e-010,        1236.334200449906,        205.1732669288739,        1739.070826217701, /* 163 HY1 */
   1.17916726160066e-009,   2.362280661261994e-010,   4.421998142341288e-010,        25.41252497179128,       -243.5820151456806,         55.6846994727801, /* 163 OXY */
  1.098855397034494e-009,   2.137550577408458e-010,   3.870185827870248e-010,        1183.074678548037,        120.6533004005955,       -1799.808420986776, /* 163 HY2 */
  2.492595392462038e-009,    3.25358042396856e-009,   3.418730889301742e-009,        185.5769608118022,        496.0838158723222,       -496.6577194284736, /* 164 HY1 */
  2.439459973186186e-009,   3.225285324428384e-009,   3.338880939017526e-009,       -142.6779273636358,       -191.8558763861422,       -36.93053477775166, /* 164 OXY */
  2.357497063885022e-009,   3.176526161622868e-009,   3.368957932573992e-009,        509.0158339063019,       -847.1086458534328,        688.2122018323811, /* 164 HY2 */
  2.861723286942781e-009,   4.972552392732735e-011,    2.08958952720236e-009,       -715.7051478963464,        435.0960706088868,        -123.759111400808, /* 165 HY1 */
  2.837479351003203e-009,   5.266798441958693e-011,   2.186561551500149e-009,       -439.3263419696864,       -335.1181297070204,       -29.54165013398794, /* 165 OXY */
  2.892843571410716e-009,  -1.375366699023441e-011,   2.236790699638084e-009,        982.6851891257464,        531.6608594210616,       -436.0982215203738, /* 165 HY2 */
  1.782690590765468e-009,    2.76505495561977e-009,   1.786901127949806e-010,       -1831.613789613948,        355.6122663324724,       -2130.259961747438, /* 166 HY1 */
  1.792210446728425e-009,   2.862676956527046e-009,   1.981661822471769e-010,        -22.6113950620343,       -292.4414177320643,        360.8690610861929, /* 166 OXY */
  1.703763243574557e-009,    2.89980831501624e-009,   2.264213494030636e-010,        627.4235185979323,        662.6182199726312,        1158.020313467634, /* 166 HY2 */
  3.479533660572121e-009,   2.353477794998776e-009,    2.50651495992505e-009,        1091.046293340362,        176.2736209877052,       -340.0602051977426, /* 167 HY1 */
  3.466209768789644e-009,   2.322361720754594e-009,   2.412417862034614e-009,       -405.7754834083815,        23.74854620703857,       -83.86602105273879, /* 167 OXY */
  3.491686362218682e-009,   2.395529707027022e-009,   2.349192954188968e-009,        867.0183033059244,       -895.6128867163603,       -645.9334626850172, /* 167 HY2 */
  2.289692273976172e-009,   1.220265129444152e-009,   1.041398031367638e-009,        1809.914314963914,        556.6941589059068,        137.0887578597089, /* 168 HY1 */
   2.27367070739279e-009,   1.196203811630611e-009,   1.137128707690823e-009,        633.1082772878334,        562.9718781702858,       -54.57200719025958, /* 168 OXY */
  2.301055528371372e-009,   1.101291129640542e-009,   1.152673944177403e-009,       -569.8830269367695,        227.4015581478224,         41.0235757180778, /* 168 HY2 */
  2.957949030181684e-009,   1.734728333501421e-009,   2.430013744769954e-009,       -2184.162212555533,        1383.838668546091,        1342.306690065802, /* 169 HY1 */
  2.901043197352666e-009,   1.711579038822573e-009,    2.35110977563082e-009,        -68.2467004520594,        270.8480745480123,        119.9922359903286, /* 169 OXY */
  2.830197884341922e-009,   1.646740757489673e-009,   2.378983399331478e-009,       -1984.367218191863,        2104.675285249992,       -418.6819214204125, /* 169 HY2 */
  7.152234189800899e-011,   1.270803590919791e-009,  -7.556577861633335e-010,        701.4704293172491,        190.2300024060278,       -737.7463094067722, /* 170 HY1 */
  1.205459662771326e-010,   1.188161159994965e-009,  -7.833509778098459e-010,        809.6732384670007,        166.2547511699847,       -475.3862267866029, /* 170 OXY */
   2.17127453620732e-010,   1.196949043918036e-009,  -7.589626732633358e-010,        799.0405355489382,         371.938433369993,       -506.9484171575251, /* 170 HY2 */
    2.5448288287924e-009,   2.921893015856743e-009,   1.160262147184314e-009,       -714.8441344618481,        888.2768494507925,       -1845.221919372269, /* 171 HY1 */
  2.507375625770362e-009,   2.836651284070044e-009,   1.123777838976605e-009,       -319.1375235099092,        21.76389727070861,       -250.5639186820081, /* 171 OXY */
  2.568440847339382e-009,   2.760874742719513e-009,   1.146776849730941e-009,         122.173673223459,        756.7370972907372,        1024.944274631742, /* 171 HY2 */
  3.943800521357092e-009,   2.044627884388196e-010,   1.653354931414908e-009,       -1377.503808520567,        865.6186168937843,        1383.143305023896, /* 172 HY1 */
  3.997587939403903e-009,   1.547727333286876e-010,   1.721456415121286e-009,       -12.62368512640817,       -176.7413849225968,       -432.4827208518602, /* 172 OXY */
   4.08945891883913e-009,   1.939521104338302e-010,   1.726426279029446e-009,        165.2920577506556,       -404.2450072483845,       -1825.987035112032, /* 172 HY2 */
  4.381018452853808e-009,   2.147481069079302e-009,   2.086798317820714e-009,        72.41869856490288,       -1766.542325833337,       -165.9016158799004, /* 173 HY1 */
  4.334333143457276e-009,   2.163389913312296e-009,   1.999807567957745e-009,       -416.8875858095712,        329.7490287841446,        465.6002944195771, /* 173 OXY */
  4.395062334900179e-009,   2.138521359263847e-009,    1.92435214849544e-009,       -1120.963807000157,         1397.78562353368,       -461.3366346777055, /* 173 HY2 */
 -6.171261675215764e-010,    3.71715434172704e-009,   2.578739835543558e-009,       -187.8956375904102,        268.1310607563414,       -336.9518878841147, /* 174 HY1 */
 -7.131030292377109e-010,   3.689296347716731e-009,   2.575222139586869e-009,       -151.3080554101406,        80.37088855954171,        133.4132284679084, /* 174 OXY */
 -7.623244778302267e-010,    3.72947709902599e-009,   2.652441024729546e-009,        129.5724060361978,        54.77477765270047,        326.1509031262886, /* 174 HY2 */
 -2.213329495449856e-010,   1.801039864990958e-009,    4.06394778019666e-010,       -989.8857500659948,        466.7057754659047,        1051.395542645496, /* 175 HY1 */
 -2.825585216774958e-010,     1.8374496710912e-009,   3.362110524360924e-010,        272.1112984100328,       -363.7446445660964,       -497.0117495915992, /* 175 OXY */
 -3.383460155425445e-010,   1.763305915132674e-009,   2.989221937699668e-010,        2519.289011440532,       -1206.566982386046,       -2242.193984327198, /* 175 HY2 */
  9.976870284426548e-010,   1.703651447623125e-009,   5.178988873943779e-009,       -2088.478853506776,       -1079.070187517661,        -603.689183697511, /* 176 HY1 */
  1.039438748721737e-009,   1.636213157096814e-009,   5.118088374564658e-009,         82.4981898118384,       -17.82665510599308,        -314.810287423027, /* 176 OXY */
  1.138728765786233e-009,   1.637589975746189e-009,   5.129903486994885e-009,       -144.9051853622766,        1417.074759125915,         1547.00331509188, /* 176 HY2 */
  2.335868319529431e-009,   2.566306168479764e-009,   3.070956547479356e-009,        691.7002220877051,        452.5484594624664,       -617.5403175481702, /* 177 HY1 */
  2.338320033737592e-009,    2.52282117396711e-009,   3.160973458944235e-009,        237.2167861005867,       -266.8825551040416,       -950.3828873183052, /* 177 OXY */
   2.37025125749512e-009,   2.588794853428738e-009,   3.229001913047214e-009,        -1633.93717133833,       -335.4717208681018,        10.69604322465537, /* 177 HY2 */
  3.148669749864374e-009,   3.476394027684496e-009,   1.858725453685536e-009,        3857.755121972385,       -1024.272458993989,       -144.8080262760082, /* 178 HY1 */
  3.157742940202709e-009,   3.570078788941588e-009,   1.892501817096742e-009,        368.1910548662862,       -602.5232996580343,       -285.6254509156498, /* 178 OXY */
  3.190399054570848e-009,   3.568422876545728e-009,   1.987004918360734e-009,       -1133.333159088085,        739.6346080837561,         268.294358719539, /* 178 HY2 */
  1.864860220045896e-009,   4.659863126257609e-009,   3.044128906012668e-009,       -1758.729121332492,        2862.481864605321,       -517.5424625685021, /* 179 HY1 */
  1.906012112691383e-009,   4.664022790856659e-009,    2.95308375829947e-009,       -138.4569369825828,       -156.5511414888386,          43.776544908165, /* 179 OXY */
  1.942349812751666e-009,   4.755804214949787e-009,    2.93709183078068e-009,       -858.1325194094201,       -339.3268901269448,       -2773.170063272119, /* 179 HY2 */
 -6.137193614678078e-010,   3.325571124015016e-010,   5.646907390049054e-010,         222.654491099011,       -2382.169811315134,        1891.259719355929, /* 180 HY1 */
 -5.558122389078582e-010,   3.981736005971184e-010,   6.130771176016326e-010,        464.4106010294496,       -1103.786453120004,       -102.3920209779363, /* 180 OXY */
 -4.642190121584741e-010,   3.977319302587508e-010,   5.729462598865259e-010,        146.4508834203639,       -1833.495223232307,       -827.2812173099752, /* 180 HY2 */
  5.119301187937802e-009,   1.893239792437662e-009,   1.863920901636871e-009,       -1161.016195117267,        307.8468012724182,        218.6916879459795, /* 181 HY1 */
  5.161607053629931e-009,   1.983830922510718e-009,   1.865781245725155e-009,       -503.8238865740819,        19.41651934978323,       -535.4256238479027, /* 181 OXY */
  5.130990451224657e-009,   2.036549762591993e-009,   1.786513676331087e-009,       -1141.699461975241,         38.6961790119698,        -277.721656710108, /* 181 HY2 */
   1.94534938743219e-010,  -1.940310267978886e-010,   2.744460735437498e-009,         1670.61320828058,       -343.3266547626102,       -992.8972284913506, /* 182 HY1 */
  1.606859031555494e-010,  -2.858348950561269e-010,   2.723813881861014e-009,       -295.4558601339969,        196.6514574480126,        -228.036930260021, /* 182 OXY */
  2.366800499074574e-010,   -3.50776551128751e-010,   2.726547188648496e-009,       -1342.587798169365,       -1103.038288721535,       -1573.956699271719, /* 182 HY2 */
   1.92293021743001e-009,   3.378845803134616e-009,   5.995585690796803e-010,        475.2167255019291,        229.8759791097748,       -1027.435623295268, /* 183 HY1 */
   1.96705418874154e-009,    3.29271024960644e-009,   6.247327934787491e-010,       -465.9457343317904,        17.44577514005209,       -86.63050832608974, /* 183 OXY */
  1.917880022457324e-009,   3.251193405575467e-009,   7.012720189253358e-010,       -1951.741609965898,        533.8205860577275,         -751.58916297622, /* 183 HY2 */
   1.89175959143061e-010,   1.865220745723826e-009,  -7.083751953342134e-011,        1414.781429439538,        127.8343097116278,        859.0760837587518, /* 184 HY1 */
  1.538106105014278e-010,   1.784734736598472e-009,  -2.317942000637428e-011,        317.4102090034114,       -423.5953676882682,       -863.0403010681531, /* 184 OXY */
  1.715049428828204e-010,   1.793589564334155e-009,   7.484354989021845e-011,        37.62046604682107,       -2555.386495526062,       -608.0051741192864, /* 184 HY2 */
  2.560832753485847e-009,   2.036849278509859e-010,    2.93350551044206e-009,       -799.2344995799672,       -1079.976617011386,        2804.461358220273, /* 185 HY1 */
  2.472896430645432e-009,   2.492464327808572e-010,   2.919672492194588e-009,        153.9136335199738,         22.8471261731796,        218.3916420778917, /* 185 OXY */
  2.488476848204156e-009,   3.435321032084262e-010,   2.890219768369437e-009,        2200.421766920304,       -172.9373923842325,        636.8712254085461, /* 185 HY2 */
 -8.563143307003503e-010,  -6.974077249877707e-010,   5.206491012828621e-010,        479.5321949443909,        1213.086099476921,         59.3983183003655, /* 186 HY1 */
 -7.882077226005627e-010,  -7.390874627544325e-010,   4.604470299493532e-010,         232.606454320034,          81.312259872434,         557.012196966958, /* 186 OXY */
 -8.327024469014004e-010,  -8.063320611336252e-010,   4.012998956298098e-010,       -306.6507057347193,        210.2702789030778,        814.4876963225335, /* 186 HY2 */
  8.634530696504649e-010,   2.681475952108643e-009,   2.802908725981249e-009,          733.47746393998,        12.36234489041196,        1010.170085973753, /* 187 HY1 */
   9.26666264765676e-010,    2.70144826292774e-009,   2.728040850601524e-009,       -128.3718198493111,       -553.5339488218688,        125.3582495659737, /* 187 OXY */
  9.197600547794758e-010,   2.798139781428692e-009,   2.703483682328108e-009,        1774.485607961938,       -46.45734132546139,        1527.295894108354, /* 187 HY2 */
  2.716821519027348e-009,   1.860928334837361e-009,   1.214904974636259e-010,        287.7496470196081,       -85.60966627422306,       -745.6806248553993, /* 188 HY1 */
  2.622285582415302e-009,   1.835422853794451e-009,   1.417987993535978e-010,         406.412890803854,        164.4995831720207,         131.226672043613, /* 188 OXY */
  2.564042779317068e-009,   1.858572241941625e-009,   6.387651938390076e-011,       -211.3983814544276,        69.64440043214472,        562.9767024754658, /* 188 HY2 */
  2.529566056833676e-010,   2.386801234512242e-009,   2.006866037763093e-009,         497.333234997233,        552.3766018562411,        610.7970038578408, /* 189 HY1 */
  2.306994124225904e-010,   2.348467887334144e-009,   2.096505157807178e-009,       -355.0443778242003,        -672.124334555959,        -116.808452148052, /* 189 OXY */
  1.378136310632248e-010,    2.31146345660684e-009,   2.094801124999455e-009,        329.8644281052188,       -2349.727076413602,       -2050.682896861507, /* 189 HY2 */
  3.755556052159544e-009,   2.065758934311146e-009,  -1.598366121673614e-009,        295.4977579367741,        942.8587437816078,        1707.641189174101, /* 190 HY1 */
  3.738143723717946e-009,   1.990097248731469e-009,  -1.535341363415576e-009,       -821.1911877376091,        106.5649626972704,        409.5529391198127, /* 190 OXY */
  3.653822228718907e-009,   2.007746400722877e-009,  -1.484563485209497e-009,       -1868.045483317488,       -698.1955713139323,       -1030.341188368004, /* 190 HY2 */
  3.491851627498078e-009,    3.02345609915244e-009,  -2.209796974721441e-009,        -299.998912894102,        1023.573979276458,       -7.083949148120557, /* 191 HY1 */
  3.419348697146465e-009,    2.97301452069051e-009,  -2.162903774226076e-009,       -296.8202620665226,        664.1285362206389,       -387.3540036445536, /* 191 OXY */
  3.453623909115166e-009,   2.939121565838554e-009,  -2.075288265005934e-009,       -299.2459050161743,        16.95288199943744,       -635.3861780203017, /* 191 HY2 */
 -1.069628146921371e-009,    2.52300770668957e-009,   3.837720912149086e-009,       -1847.072013632158,        -1514.75350512221,        335.5322616612691, /* 192 HY1 */
 -1.067124091698629e-009,   2.535017322305736e-009,   3.936965553516736e-009,        -737.713746394128,        58.73683929807658,        126.5803610041023, /* 192 OXY */
 -1.152386967115329e-009,   2.501195875252192e-009,   3.976794585577296e-009,       -514.7267482775417,        1169.470855119696,        1568.232208977347, /* 192 HY2 */
   1.84085013141135e-010,   3.036183876411686e-009,    2.59917116766967e-009,         124.098197493322,       -426.6185975640808,        85.05531382461318, /* 193 HY1 */
  2.499475228072122e-010,   2.964821963479424e-009,   2.575304722539849e-009,        128.4846123588164,       -320.6787840800523,       -220.7024135228183, /* 193 OXY */
  2.016742776259918e-010,   2.878498145912658e-009,   2.560543546422244e-009,        122.6792973261912,        -237.654854436466,       -691.1089275927063, /* 193 HY2 */
  1.092059837133282e-009,   2.672784669675992e-009,   3.532135323993844e-009,       -1730.858286639971,        714.9245285633083,         933.904921083495, /* 194 HY1 */
  1.126422824639487e-009,   2.640324899716574e-009,   3.620257673167198e-009,        443.3118812817972,       -431.5307003239448,       -314.6386910821812, /* 194 OXY */
  1.120966493190143e-009,    2.71420899892616e-009,   3.687424393580694e-009,       -53.93209714784969,       -1966.411565347158,         1353.39758378356, /* 194 HY2 */
  1.297067284288506e-010,   2.428784041252907e-009,  -1.871056066629353e-009,        -206.922336821805,       -542.6479556017098,       -861.4985910045302, /* 195 HY1 */
  2.194038653150779e-010,   2.470430448476742e-009,  -1.856223656179238e-009,       -545.1869126412046,       -373.6512687238281,        756.1196009378692, /* 195 OXY */
  2.124181227842927e-010,   2.569541816295462e-009,   -1.86754338101677e-009,       -70.08840716696953,        -601.683631495993,       -1669.747923186446, /* 195 HY2 */
  1.670379239802431e-009,   2.137909198856423e-009,   2.389829379856766e-009,       -2921.954781783911,        -546.380241221549,        681.3742853838473, /* 196 HY1 */
  1.679601078024821e-009,   2.040336913344177e-009,   2.409694092281494e-009,        39.44344257590622,       -351.9493786216188,        373.6488508600084, /* 196 OXY */
  1.710741140218058e-009,   2.028310540366177e-009,   2.503957886389284e-009,         3383.20825661395,          334.73398467823,        -609.887032660145, /* 196 HY2 */
  1.324296462033267e-009,    3.26085198686053e-009,    -3.4839766065727e-011,        1084.449232856086,         -102.42528248381,       -222.6071276047323, /* 197 HY1 */
  1.350471222537721e-009,   3.338005540546684e-009,  -9.282434217121654e-011,        431.3215267122168,        205.5457172063802,       -109.9553136615827, /* 197 OXY */
  1.304658308507748e-009,   3.330202971966372e-009,  -1.813697928258509e-010,       -538.0935191059659,        656.5854584527186,        348.0492106728294, /* 197 HY2 */
 -3.732604788042243e-010,   2.651689846166314e-010,   3.411863913653648e-009,       -1597.427792830975,       -1105.738719912831,       -803.8393372876891, /* 198 HY1 */
 -4.387548785700859e-010,   2.677430322099482e-010,   3.487387803107892e-009,        -300.028978046186,        472.7659809392318,        285.2127468084954, /* 198 OXY */
 -5.313817039416926e-010,   2.770824461945332e-010,   3.450876892848048e-009,       -877.6496855414755,         393.391637017668,        1714.005228142557, /* 198 HY2 */
  3.034788499821506e-009,    3.41708879647584e-009,   2.296034154591485e-009,        770.0637688661162,        1248.438450652094,        799.7455345334575, /* 199 HY1 */
  3.119924166123926e-009,   3.369606147785039e-009,   2.273733562748458e-009,       -239.1158459194364,        -41.0388525301041,       -352.3085121556668, /* 199 OXY */
  3.162036166917156e-009,   3.335452703864564e-009,   2.357758092842726e-009,         604.919437720437,       -289.8981046998552,        -873.370124464354, /* 199 HY2 */
 -5.107800595182579e-010,   1.424129418708969e-009,    2.62748987151904e-009,        -2592.74015050747,       -444.6968297739502,       -17.68430329977981, /* 200 HY1 */
 -5.227111163147348e-010,   1.472616274030794e-009,   2.540848803054248e-009,        316.9520056002228,        321.9460699781138,       -15.45802298436894, /* 200 OXY */
 -6.074214491949376e-010,   1.525680410844676e-009,   2.543739631226976e-009,        461.3864680899567,        695.1972758422917,       -2205.753693712186, /* 200 HY2 */
  2.267827218785277e-009,   2.411627040777922e-009,   3.362380639701145e-009,        281.3822930578458,       -520.5037307380753,        311.4107660542882, /* 201 HY1 */
  2.237386528184472e-009,   2.318946524119226e-009,   3.384373506445235e-009,       -128.0792225319151,        -427.581161511446,        138.5983600097248, /* 201 OXY */
  2.137690833941244e-009,   2.317547450498613e-009,   3.392042348754534e-009,       -58.60861343633643,        230.2962714583674,        1213.689164136407, /* 201 HY2 */
   1.88999119990326e-009,   1.292217082427654e-009,    3.26169525753575e-010,         287.099699116063,       -1631.807770439836,       -604.3317290344719, /* 202 HY1 */
  1.799332551671723e-009,     1.3154598997536e-009,    2.90944805121846e-010,        542.0268932256582,        48.39578439779731,       -173.5857721346955, /* 202 OXY */
  1.730024454972069e-009,   1.262278871767728e-009,   3.396086109146282e-010,       -146.2769676979608,        721.4698671577991,       -413.2751726286764, /* 202 HY2 */
   1.23707565667276e-009,   1.083612507860757e-009,   3.207817827111688e-009,       -1219.430950160311,        63.45995031602154,       -879.5725288719216, /* 203 HY1 */
  1.302832069937195e-009,   1.105514134734565e-009,   3.135731678269642e-009,       -401.0698519480381,        313.5680205757703,       -61.93881287428322, /* 203 OXY */
  1.394818494739504e-009,   1.080234266987079e-009,   3.165722109056126e-009,       -748.3477863640892,        134.5003841681804,        860.6558525713939, /* 203 HY2 */
  1.137841033894716e-009,  -4.002827440204698e-010,   4.876845265118735e-009,        61.23145012784722,        1024.958508068025,        583.9192771290345, /* 204 HY1 */
  1.042572923968121e-009,  -3.708796925977637e-010,   4.884555504936556e-009,       -218.8513673137024,         403.928605560999,       -458.2319096986357, /* 204 OXY */
  1.030074196789696e-009,    -3.2041684132213e-010,   4.969979637709765e-009,       -769.3345823902405,        3016.559238962853,       -2053.828220043712, /* 204 HY2 */
  1.029814484154682e-010,   3.477725675323058e-009,   4.014519228699734e-010,         942.096405835792,        705.6609303993037,       -683.9698860875574, /* 205 HY1 */
  1.971641897124474e-010,   3.500615286424114e-009,   3.768414180919283e-010,        1257.514988885016,        213.4546906163847,        56.28891086877613, /* 205 OXY */
  2.124681233184218e-010,   3.598425810473898e-009,   3.909442901856504e-010,        1810.617122222816,        41.28838293395978,        662.6075942953428, /* 205 HY2 */
   2.30996858345907e-009,   3.335445683140838e-009,   1.328104776724508e-009,        445.8027804723364,        1294.505669571905,       -918.0006968519132, /* 206 HY1 */
   2.25315489453733e-009,   3.311495495058794e-009,   1.249373631413778e-009,        115.5307972200662,       -198.6580405558358,       -234.3573334241587, /* 206 OXY */
  2.276646703408992e-009,   3.219468971342024e-009,   1.218080599465016e-009,         1957.30876137238,        293.8667467440208,       -329.2532546363222, /* 206 HY2 */
 -5.816052165003497e-010,   5.532896292047869e-010,     2.6799675524199e-011,        1507.046318245996,        86.38642665729547,        1814.665617422165, /* 207 HY1 */
   -5.0478631747874e-010,   5.203349892002095e-010,  -2.808974472484918e-011,        398.2493740372537,        154.6755040365979,        204.4514267226004, /* 207 OXY */
 -5.393996351916318e-010,   4.665942919757775e-010,  -1.049913422293121e-010,       -1149.643751937186,        -196.116292180423,        1135.292156534992, /* 207 HY2 */
   1.66589201474724e-009,   1.437112814804868e-009,  -1.803209041855718e-009,        35.73140123417154,       -3052.621547290132,        266.9547592199466, /* 208 HY1 */
  1.619978911779109e-009,   1.427339525315879e-009,  -1.714911421341888e-009,        189.0712464420377,       -144.9484043879209,        693.0449132499951, /* 208 OXY */
  1.679865123106499e-009,   1.378520463688997e-009,  -1.651426455289076e-009,       -1104.297665444075,       -1247.182723921909,        1077.447592977803, /* 208 HY2 */
  9.432625080066986e-010,   1.682753507175223e-009,   1.147366169247655e-009,        59.20089523095329,       -369.0397439962174,       -1078.767644399989, /* 209 HY1 */
  8.958702837441049e-010,   1.765055013404145e-009,   1.116054204251873e-009,        169.6831831227585,       -161.0756121099918,       -700.9501054801076, /* 209 OXY */
  8.918843587188097e-010,   1.765243686127503e-009,   1.016133851949142e-009,       -2654.964427937891,       -1293.933117070114,       -613.6036672940486, /* 209 HY2 */
  5.443794213302345e-012,   1.720753351521082e-009,   1.438665007165652e-009,       -1193.355532454212,       -407.0668414079134,        -312.822343344071, /* 210 HY1 */
  6.568397451354104e-011,   1.723750474520152e-009,   1.358901996119659e-009,       -569.3349223278584,       -334.6754348041588,        159.2474355261676, /* 210 OXY */
  9.527480287996592e-011,   1.631020402192354e-009,   1.335977848094701e-009,         1040.60689166038,       -134.9224454437536,        1384.298916521998, /* 210 HY2 */
  1.791709436532572e-009,   4.639056600422866e-009,   3.505374562259587e-009,        566.9815352880776,         993.251013137241,         553.613768929467, /* 211 HY1 */
  1.844285005736869e-009,   4.587424427685124e-009,    3.43777331773536e-009,       -173.5136419059749,        512.0348598521316,        342.1994745888326, /* 211 OXY */
  1.823213238516478e-009,    4.62141481637507e-009,    3.34611834458651e-009,       -253.0367093003382,         684.595160541631,        424.3594511800035, /* 211 HY2 */
  2.805257748436973e-009,   3.101086236297332e-009,   4.737810771745994e-011,        2317.363952638586,         1910.55688552278,       -1304.957019724571, /* 212 HY1 */
  2.844558266633818e-009,   3.039416383971852e-009,   -2.08296527179532e-011,       -258.8348363754422,       -741.8730575403677,       -443.9751244955687, /* 212 OXY */
  2.790089286150347e-009,   3.043311241466588e-009,  -1.046029196363904e-010,       -960.9240806737934,       -770.5479382001673,         9.10071758539093, /* 212 HY2 */
  1.380375984101732e-009,   3.430836854255707e-009,   3.046567063145678e-010,         397.482814596249,         1402.51247344563,        1567.606636600424, /* 213 HY1 */
  1.442285637007052e-009,    3.47772553231576e-009,     2.4165950997278e-010,        119.0713861003371,       -105.2438197519496,        154.5284766285749, /* 213 OXY */
  1.500619881756572e-009,   3.410515536527714e-009,    1.96053234574771e-010,       -2763.530813479924,       -1573.908587948381,        -1439.49090033326, /* 213 HY2 */
  3.346916881275488e-009,   2.279766002630146e-009,   2.845378625500042e-009,       -1454.338631151462,       -1259.478507000013,         6.61855884387563, /* 214 HY1 */
  3.274530089740992e-009,   2.319635036355606e-009,   2.789070432835562e-009,        -310.355761599667,       -68.76090881200506,          -634.8687720142, /* 214 OXY */
  3.281683054263768e-009,   2.284423062909547e-009,   2.695748666754062e-009,       -2116.972771744439,       -2840.984496450603,        241.2785483363914, /* 214 HY2 */
  1.807616970254674e-009,   5.215339720532661e-009,   2.257551142642602e-009,        304.3389086939869,        1140.112018503638,       -1993.243981362014, /* 215 HY1 */
   1.78623449835369e-009,   5.304444780073802e-009,   2.297589602397862e-009,       -6.830571882096374,        4.302082082801322,        413.1226130882086, /* 215 OXY */
  1.725428971973262e-009,   5.292026369961767e-009,   2.376001789281601e-009,       -104.7865309042549,       -2307.541761863651,       -11.32792254562439, /* 215 HY2 */
};

}
