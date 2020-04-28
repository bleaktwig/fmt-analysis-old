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
    String infile;
    boolean[] pltLArr;
    boolean debugInfo;
    boolean testRun;

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

        // Set geometry parameters reading from database
        DatabaseConstantProvider dbProvider =
                new DatabaseConstantProvider(10, "rgf_spring2020");
        String fmtTable = "/geometry/fmt/fmt_layer_noshim";
        dbProvider.loadTable(fmtTable);
        int ln = 3;

        double[] fmtZ     = new double[ln]; // z position of the layers in mm
        double[] fmtAngle = new double[ln]; // strip angle in deg

        for (int li=0; li<ln; li++) {
            fmtZ[li]     = dbProvider.getDouble(fmtTable+"/Z",li);
            fmtAngle[li] = dbProvider.getDouble(fmtTable+"/Angle",li);
        }

        DataGroup[] dgFMT = Data.createDataGroups(ln, zn, r);

        double[] absMeanArr    = new double[zn];
        double[] absMeanErrArr = new double[zn];
        double[] stddevArr     = new double[zn];
        double[] stddevErrArr  = new double[zn];

        // === RUN =============================================================
        for (int zi=0; zi<zn; ++zi) {
            int[] stopCount = new int[]{0, 0, 0, 0, 0, 0};
            int nevent = 0;
            HipoDataSource reader = new HipoDataSource();
            reader.open(infile);
            System.out.printf("\n\n");

            // === LOOP THROUGH EVENTS =========================================
            while (reader.hasEvent()) {
                if (nevent == 10000 && testRun) break;
                DataEvent event = reader.getNextEvent();
                if (nevent%50000 == 0) {
                    if (nevent == 0) {
                        System.out.printf("z shift : %d/%d\n", zi+1, zn);
                        System.out.format("Analyzed %8d events...\n", nevent);
                    }
                    else {
                        System.out.format("Analyzed %8d events...\n", nevent);
                    }
                }
                nevent++;

                // Get relevant data banks.
                // DataBank hits      = getBank(event, "FMTRec::Hits");
                DataBank clusters = getBank(event, "FMTRec::Clusters");
                DataBank traj     = getBank(event, "REC::Traj");
                DataBank particle = getBank(event, "REC::Particle");
                // DataBank recTrack = null;

                // if (event.hasBank("REC::Track"))
                //     recTrack = event.getBank("REC::Track");

                // Ignore events that don't have the necessary banks.
                if (clusters==null || traj==null || particle==null) continue;

                // Loop through trajectory points.
                for (int loop=0; loop<traj.rows(); loop++) {
                    int detector = traj.getByte("detector", loop);
                    int layer    = traj.getByte("layer", loop);
                    int pindex   = traj.getShort("pindex", loop);
                    // int sector;

                    // Use only FMT layers 1,2,3 (ignore 4,5,6 that are not
                    // installed in RG-F).
                    if (detector!=DetectorType.FMT.getDetectorId() ||
                            layer<1 || layer>ln) {
                        stopCount[0]++;
                        continue;
                    }

                    // for (int j = 0; j<recTrack.rows(); ++j) {
                    //     if (recTrack.getShort("pindex", j) == pindex) {
                    //         sector = recTrack.getByte("sector", j);
                    //     }
                    // }

                    // TODO: Make 6 different plots, one for each sector!

                    // === TRANSLATE TRAJECTORY USING SWIMMER ==================
                    double zShift = zShArr[zi]; // shift in z.
                    double zRef = fmtZ[layer-1]/10 + zShift;
                    double phiRef = fmtAngle[layer-1];

                    // Get relevant data.
                    double x  = (double) particle.getFloat("vx", pindex);
                    double y  = (double) particle.getFloat("vy", pindex);
                    double z  = (double) particle.getFloat("vz", pindex);
                    double px = (double) particle.getFloat("px", pindex);
                    double py = (double) particle.getFloat("py", pindex);
                    double pz = (double) particle.getFloat("pz", pindex);
                    int q     = (int)    particle.getByte("charge", pindex);

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
                        if (layer!=clusters.getByte("layer", i)) continue;
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
                        dgFMT[zi].getH1F("hi_cluster_res_l"+layer)
                                .fill(yLoc - yclus);
                        dgFMT[zi].getH2F("hi_cluster_res_strip_l"+layer)
                                .fill(yLoc - yclus, strip);
                    }
                }
            }
            System.out.format("Analyzed %8d events... Done!\n", nevent);
            reader.close();

            if (debugInfo) printDebugInfo(stopCount);

            // Fit residual plots
            for (int li=1; li<=ln; ++li) {
                Data.fitRes(dgFMT[zi].getH1F("hi_cluster_res_l"+li),
                        dgFMT[zi].getF1D("f1_res_l"+li), g);
            }

            // Get absolute average and standard deviation
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

        EmbeddedCanvasTabbed fmtCanvas = Data.drawPlots(dgFMT, zn, titleArr, pltLArr);
        JFrame frame = new JFrame("FMT");
        frame.setSize(1600, 1000);
        frame.add(fmtCanvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);

        return 0;
    }
}
