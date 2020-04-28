/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import javax.swing.JFrame;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.*;
import org.jlab.groot.math.F1D;

public class ResolutionAnalysis {
    // Constants:
    final int ln = 3;  // Number of FMT layers.

    // Class variables:
    String infile;
    boolean[] pltLArr;
    boolean debugInfo;
    boolean testRun;
    double[] fmtZ;     // z position of the layers in mm
    double[] fmtAngle; // strip angle in deg

    /**
     *
     * @param infile  : Input hipo file.
     * @param pltLArr : Boolean array describing what lines should be drawn in
     *                  each plot:
     *                  - [0] : Top plots, vertical line at 0.
     *                  - [1] : Bottom plots, vertical line at 0.
     *                  - [2] : Bottom plots, horizontal lines at each cable's
     *                          endpoint.
     *                  - [3] : Bottom plots, horizonal lines separating each
     *                          "region" of the FMT.
     * @param dbgInfo : Boolean describing if debugging info should be printed.
     * @param testRun : Boolean describing if the run should be cut short for
     *                  expedient testing.
     * @return status int.
     */
    public ResolutionAnalysis(String infile, boolean[] pltLArr,
            boolean debugInfo, boolean testRun) {

        // Sanitize input.
        if (pltLArr.length != 4) {
            System.out.printf("pltLArr should have a size of 4. Read the ");
            System.out.printf("method's description.\n");
            System.exit(1);
        }

        this.infile    = infile;
        this.pltLArr   = pltLArr;
        this.debugInfo = debugInfo;
        this.testRun   = testRun;

        // Set geometry parameters by reading from database
        DatabaseConstantProvider dbProvider =
                new DatabaseConstantProvider(10, "rgf_spring2020");
        String fmtTable = "/geometry/fmt/fmt_layer_noshim";
        dbProvider.loadTable(fmtTable);

        fmtZ     = new double[ln]; // z position of the layers in mm
        fmtAngle = new double[ln]; // strip angle in deg

        for (int li=0; li<ln; li++) {
            fmtZ[li]     = dbProvider.getDouble(fmtTable+"/Z",li);
            fmtAngle[li] = dbProvider.getDouble(fmtTable+"/Angle",li);
        }

        return;
    }

    private int printDebugInfo(int[] sc) {
        System.out.printf("\n");
        System.out.printf("        wrong detector or layer id : %8d\n", sc[0]);
        System.out.printf("tracks further downstream than FMT : %8d\n", sc[1]);
        System.out.printf("       tracks not in the FMT layer : %8d\n", sc[2]);
        System.out.printf("         tracks too close to (0,0) : %8d\n", sc[3]);
        System.out.printf("  clusters with wrong strip number : %8d\n", sc[4]);
        System.out.printf("      clusters with too low energy : %8d\n", sc[5]);
        System.out.printf("              TOTAL RELEVANT STOPS : %8d\n",
                sc[1] + sc[2] + sc[3] + sc[4] + sc[5]);
        System.out.printf("\n");

        return 0;
    }

    private DataBank getBank(DataEvent event, String name) {
        return (event.hasBank(name)) ? event.getBank(name) : null;
    }

    /**
     * Generic function for running analysis.
     * @param type  : Type of analysis to be ran:
     *                  * 0 : z shift.
     *                  * 1 : dc sector.
     * @param opt   : Variable used in different manner by different types:
     *                  * z shift   : plot index (zi).
     *                  * dc sector : number of sector (sn).
     * @param swim  : TrkSwim class instance.
     * @param dgFMT : Array of data groups where analysis data is stored.
     * @param zSh   : z shift to be used.
     * @param g     : Range for the gaussian fit.
     * @return
     */
    private int runAnalysis(int type, int opt, TrkSwim swim, DataGroup[] dgFMT,
            double zSh, int g) {
        int[] stopCount = new int[]{0, 0, 0, 0, 0, 0};
        int ei = 0; // Event number.
        HipoDataSource reader = new HipoDataSource();
        reader.open(infile);
        System.out.printf("\nRunning analysis...\n");

        // === LOOP THROUGH EVENTS =========================================
        while (reader.hasEvent()) {
            if (ei == 10000 && testRun) break;
            if (ei%50000==0) System.out.format("Analyzed %8d events...\n", ei);
            DataEvent event = reader.getNextEvent();
            ei++;

            // Get relevant data banks.
            DataBank clusters = getBank(event, "FMTRec::Clusters");
            DataBank traj     = getBank(event, "REC::Traj");
            DataBank particle = getBank(event, "REC::Particle");
            DataBank track    = null;
            if (clusters==null || traj==null || particle==null) continue;

            if (type == 1) {
                track = getBank(event, "REC::Track");
                if (track==null) continue;
            }

            // Loop through trajectory points.
            for (int loop=0; loop<traj.rows(); loop++) {
                int detector = traj.getByte("detector", loop);
                int li = traj.getByte("layer", loop);
                int pi = traj.getShort("pindex", loop);
                int si = -1; // DC sector.

                // Use only FMT layers 1, 2, and 3.
                if (detector!=DetectorType.FMT.getDetectorId() || li<1 || li>ln) {
                    stopCount[0]++;
                    continue;
                }

                if (type == 1) {
                    for (int row = 0; row<track.rows(); ++row) {
                        if (track.getShort("pindex", row) == pi) {
                            si = track.getByte("sector", row)-1;
                        }
                    }
                }

                // === TRANSLATE TRAJECTORY USING SWIMMER ==================
                double zRef = fmtZ[li-1]/10 + zSh;
                double phiRef = fmtAngle[li-1];

                // Get relevant data.
                double x  = (double) particle.getFloat("vx", pi);
                double y  = (double) particle.getFloat("vy", pi);
                double z  = (double) particle.getFloat("vz", pi);
                double px = (double) particle.getFloat("px", pi);
                double py = (double) particle.getFloat("py", pi);
                double pz = (double) particle.getFloat("pz", pi);
                int q     = (int)    particle.getByte("charge", pi);

                // See if the track is too downstream.
                if (z > zRef) {
                    stopCount[1]++;
                    continue;
                }

                double[] V = swim.swimToPlane(x,y,z,px,py,pz,q,zRef);

                x  = V[0];
                y  = V[1];
                z  = V[2];
                px = V[3];
                py = V[4];
                pz = V[5];

                // Ignore tracks that aren't in the FMT layer.
                if (Math.abs(z-zRef)>0.05) {
                    stopCount[2]++;
                    continue;
                }
                if (25.0 > x*x + y*y || x*x + y*y > 225.0) {
                    stopCount[3]++;
                    continue;
                }

                // Rotate (x,y) to local coordinates.
                double xLoc = x * Math.cos(Math.toRadians(phiRef))
                        + y * Math.sin(Math.toRadians(phiRef));
                double yLoc = y * Math.cos(Math.toRadians(phiRef))
                        - x * Math.sin(Math.toRadians(phiRef));

                // Loop over the clusters and calculate residuals for every
                // track-cluster combination.
                for (int i=0; i<clusters.rows(); i++) {
                    // Check that the cluster and trajectory layers match.
                    if (li!=clusters.getByte("layer", i)) continue;
                    int strip     = clusters.getInt("seedStrip", i);
                    double yclus  = clusters.getFloat("centroid", i);
                    double energy = clusters.getFloat("ETot", i);

                    // Check the strip number.
                    if (strip<0 || strip>1023) {
                        stopCount[4]++;
                        continue;
                    }
                    // Apply energy cuts.
                    if (energy<=50) {
                        stopCount[5]++;
                        continue;
                    }
                    int plti = -1;
                    if (type == 0) plti = opt;
                    if (type == 1) plti = si;
                    dgFMT[plti].getH1F("hi_cluster_res_l"+li)
                            .fill(yLoc - yclus);
                    dgFMT[plti].getH2F("hi_cluster_res_strip_l"+li)
                            .fill(yLoc - yclus, strip);
                }
            }
        }
        System.out.format("Analyzed %8d events... Done!\n", ei);
        reader.close();

        if (debugInfo) printDebugInfo(stopCount);

        // Fit residual plots
        if (type == 0) {
            for (int li=1; li<=ln; ++li) {
                Data.fitRes(dgFMT[opt].getH1F("hi_cluster_res_l"+li),
                        dgFMT[opt].getF1D("f1_res_l"+li), g);
            }
        }
        else if(type == 1) {
            // Fit residual plots
            for (int si=0; si<opt; ++si) {
                for (int li=1; li<=ln; ++li) {
                    Data.fitRes(dgFMT[si].getH1F("hi_cluster_res_l"+li),
                            dgFMT[si].getF1D("f1_res_l"+li), g);
                }
            }
        }

        return 0;
    }

    /**
     *
     * @param zShArr  : Array of global z shifts to try.
     * @param r       : Plot range.
     * @param g       : Gaussian range.
     * @param swim    : TrkSwim class instance.
     * @return status int.
     */
    public int zShiftAnalysis(double[] zShArr, int r, int g, TrkSwim swim) {
        // === SETUP ===========================================================
        // Set canvases' stuff.
        int zn = zShArr.length;
        String[] titleArr = new String[zn];
        for (int zi=0; zi < zn; ++zi) titleArr[zi] = "z shift : "+zShArr[zi];

        DataGroup[] dgFMT = Data.createDataGroups(ln, zn, r);

        double[] absMeanArr    = new double[zn];
        double[] absMeanErrArr = new double[zn];
        double[] stddevArr     = new double[zn];
        double[] stddevErrArr  = new double[zn];

        // === RUN =============================================================
        for (int zi=0; zi<zn; ++zi) {
            runAnalysis(0, zi, swim, dgFMT, zShArr[zi], g);

            // Get absolute mean and standard deviation from fit.
            for (int li=1; li<=ln; ++li) {
                F1D func = dgFMT[zi].getF1D("f1_res_l"+li);
                absMeanArr[zi]    += Math.abs(func.getParameter(1));
                absMeanErrArr[zi] += func.parameter(1).error();
                stddevArr[zi]     += func.getParameter(2);
                stddevErrArr[zi]  += func.parameter(2).error();
            }

            absMeanArr[zi]    /= ln;
            absMeanErrArr[zi] /= ln;
            stddevArr[zi]     /= ln;
            stddevErrArr[zi]  /= ln;
        }

        // === PRINT ALIGNMENT DATA AND DRAW PLOTS =============================
        for (int zi=0; zi<zn; ++zi) {
            System.out.printf("z shift : %5.2f\n", zShArr[zi]);
            System.out.printf("  * |mean| : %9.6f +- %9.6f\n",
                    absMeanArr[zi], absMeanErrArr[zi]);
            System.out.printf("  * omega  : %9.6f +- %9.6f\n",
                    stddevArr[zi], stddevErrArr[zi]);
        }

        Data.drawPlots(dgFMT, zn, titleArr, pltLArr);

        return 0;
    }

    /**
     *
     * @param zSh  : z shift.
     * @param r    : Plot range.
     * @param g    : Gaussian range.
     * @param swim : TrkSwim class instance.
     * @return status int.
     */
    public int dcSectorAnalysis(double zSh, int r, int g, TrkSwim swim) {
        // === SETUP ===========================================================
        // Set canvases' stuff.
        int sn = 6; // Number of DC sectors.
        String[] titleArr = new String[sn];
        for (int si=0; si < sn; ++si) titleArr[si] = "DC sector "+(si+1);

        DataGroup[] dgFMT = Data.createDataGroups(ln, sn, r);
        runAnalysis(1, sn, swim, dgFMT, zSh, g);
        Data.drawPlots(dgFMT, sn, titleArr, pltLArr);

        return 0;
    }
}
