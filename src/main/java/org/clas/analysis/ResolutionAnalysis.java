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

public class ResolutionAnalysis {
    /**
     *
     * @param infile    : Input hipo file.
     * @param zShArr    : Array of global z shifts to try.
     * @param pltRanArr : Array describing plotting ranges. One plot is produced
     *                    for every z shift times plot range.
     * @param gssRanArr : Array describing ranges for gaussian fitting.
     * @param pltLArr   : Boolean array describing what lines should be drawn
     *                    in each plot:
     *                    - [0] : Top plots, vertical line at 0.
     *                    - [1] : Bottom plots, vertical line at 0.
     *                    - [2] : Bottom plots, horizontal lines at each cable's
     *                            endpoint.
     *                    - [3] : Bottom plots, horizonal lines separating each
     *                            "region" of the FMT.
     * @param swim      : TrkSwim class instance.
     * @param dbgInfo   : Boolean describing if debugging info should be printed.
     * @return status int.
     */
    public int runAnalysis(String infile, double[] zShArr, int[] pltRanArr,
            int[] gssRanArr, boolean[] pltLArr, TrkSwim swim, boolean dbgInfo) {

        // Sanitize input.
        if (pltRanArr.length != gssRanArr.length) {
            System.out.printf("pltRanArr and gssRanArr should have the same ");
            System.out.printf("size. Read the method's description.\n");
            return 1;
        }
        if (pltLArr.length != 4) {
            System.out.printf("pltLArr should have a size of 4. Read the ");
            System.out.printf("method's description.\n");
            return 1;
        }

        // Set canvases' stuff
        int zn = zShArr.length;
        int rn = pltRanArr.length;
        int cn = zn * rn;
        double[] avgMeanArr = new double[3];
        String[] titleArr = new String[cn];
        for (int zi=0; zi < zn; ++zi) {
            for (int ri=0; ri < rn; ++ri) {
                titleArr[zi * rn + ri] = "z shift : "+zShArr[zi] +
                        ", range : [-"+pltRanArr[ri]+","+pltRanArr[ri]+"]";
            }
        }

        // === SET GEOMETRY FROM DB ============================================
        // Set geometry parameters reading from database
        DatabaseConstantProvider dbProvider =
                new DatabaseConstantProvider(10, "rgf_spring2020");
        String fmtTable = "/geometry/fmt/fmt_layer_noshim";
        dbProvider.loadTable(fmtTable);
        int ln = 3;
        // z position of the FMT layers in mm
        double[] fmtZ     = new double[ln];
        // strip angle in deg
        double[] fmtAngle = new double[ln];

        for (int li=0; li<ln; li++) {
            fmtZ[li]     = dbProvider.getDouble(fmtTable+"/Z",li);
            fmtAngle[li] = dbProvider.getDouble(fmtTable+"/Angle",li);
        }

        DataGroup[] dgFMT = Data.createDataGroups(ln, zn, rn, pltRanArr);

        // === OPEN INPUT FILE =================================================
        for (int zi=0; zi<zn; ++zi) {
            int[] stopCount = new int[]{0, 0, 0, 0, 0, 0, 0};
            int nevent = 0;
            // HipoDataSource : reader object, get data from bank.
            HipoDataSource reader = new HipoDataSource();
            reader.open(infile);
            System.out.printf("\n\n");

            // === LOOP THROUGH EVENTS =========================================
            while (reader.hasEvent()) {
                // FOR QUICKER TESTING, REMOVE BEFORE PRODUCTION!
                if (nevent == 10000) break;
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
                DataBank fmtHits = null;
                DataBank fmtClusters = null;
                DataBank recTrajEB = null;
                DataBank recParticles = null;
                if (event.hasBank("FMTRec::Hits"))
                    fmtHits = event.getBank("FMTRec::Hits");
                if (event.hasBank("FMTRec::Clusters"))
                    fmtClusters = event.getBank("FMTRec::Clusters");
                if (event.hasBank("REC::Traj"))
                    recTrajEB = event.getBank("REC::Traj");
                if (event.hasBank("REC::Particle"))
                    recParticles = event.getBank("REC::Particle");

                // Ignore events that don't have the necessary banks.
                if (fmtHits==null || fmtClusters==null ||
                        recTrajEB==null || recParticles==null) {
                    stopCount[0]++;
                    continue;
                }

                // Loop through trajectory points.
                for (int loop=0; loop<recTrajEB.rows(); loop++) {
                    int detector = recTrajEB.getByte("detector", loop);
                    int layer    = recTrajEB.getByte("layer", loop);
                    int pindex   = recTrajEB.getShort("pindex", loop);

                    // Use only FMT layers 1,2,3 (ignore 4,5,6 that are not
                    // installed in RG-F).
                    if (detector!=DetectorType.FMT.getDetectorId() ||
                            layer<1 || layer>ln) {
                        stopCount[1]++;
                        continue;
                    }

                    // === TRANSLATE TRAJECTORY USING SWIMMER ==================
                    double zShift = zShArr[zi]; // shift in z.
                    double zRef = fmtZ[layer-1]/10 + zShift;
                    double phiRef = fmtAngle[layer-1];

                    // Get relevant data.
                    double x  = (double) recParticles.getFloat("vx", pindex);
                    double y  = (double) recParticles.getFloat("vy", pindex);
                    double z  = (double) recParticles.getFloat("vz", pindex);
                    double px = (double) recParticles.getFloat("px", pindex);
                    double py = (double) recParticles.getFloat("py", pindex);
                    double pz = (double) recParticles.getFloat("pz", pindex);
                    int q     = (int)    recParticles.getByte("charge", pindex);

                    // See if the track is too downstream.
                    if (z > zRef) {
                        stopCount[2]++;
                        continue;
                    }

                    double[] V = swim.swimToPlane(x,y,z,px,py,pz,q,zRef);

                    x  = V[0];
                    y  = V[1];
                    z  = V[2];
                    px = V[3];
                    py = V[4];
                    pz = V[5];

                    // Ignore tracks that aren't in the FMT layers.
                    if (Math.abs(z-zRef)>0.05) {
                        stopCount[3]++;
                        continue;
                    }
                    // Ignore tracks that go through the hole at the center of
                    // the FMT layers.
                    if (25.0 > x*x + y*y || x*x + y*y > 225.0) {
                        stopCount[4]++;
                        continue;
                    }

                    // Rotate (x,y) to local coordinates.
                    double xLoc = x * Math.cos(Math.toRadians(phiRef))
                            + y * Math.sin(Math.toRadians(phiRef));
                    double yLoc = y * Math.cos(Math.toRadians(phiRef))
                            - x * Math.sin(Math.toRadians(phiRef));

                    // Loop over the clusters and calculate residuals for every
                    // track-cluster combination.
                    for (int i=0; i<fmtClusters.rows(); i++) {
                        // Check that the cluster layer matches the trajectory
                        // layer
                        if (layer!=fmtClusters.getByte("layer", i)) continue;
                        int strip     = fmtClusters.getInt("seedStrip", i);
                        double yclus  = fmtClusters.getFloat("centroid", i);
                        double energy = fmtClusters.getFloat("ETot", i);

                        // apply minimal cuts (to be optimized)
                        if (strip<0 || strip>1023) {
                            stopCount[5]++;
                            continue;
                        }
                        if (energy<=50) {
                            stopCount[6]++;
                            continue;
                        }
                        for (int ri=0; ri<rn; ++ri) {
                            dgFMT[zi*rn + ri]
                                    .getH1F("hi_cluster_res_l" + layer)
                                    .fill(yLoc - yclus);
                            dgFMT[zi*rn + ri]
                                    .getH2F("hi_cluster_res_strip_l" + layer)
                                    .fill(yLoc - yclus, strip);
                        }
                    }
                }
            }
            System.out.format("Analyzed %8d events... Done!\n", nevent);
            reader.close();

            if (dbgInfo) {
                System.out.printf("\n");
                System.out.printf("    events without necessary files : %8d\n",
                        stopCount[0]);
                System.out.printf("        wrong detector or layer id : %8d\n",
                        stopCount[1]);
                System.out.printf("tracks further downstream than FMT : %8d\n",
                        stopCount[2]);
                System.out.printf("       tracks not in the FMT layer : %8d\n",
                        stopCount[3]);
                System.out.printf("         tracks too close to (0,0) : %8d\n",
                        stopCount[4]);
                System.out.printf("  clusters with wrong strip number : %8d\n",
                        stopCount[5]);
                System.out.printf("      clusters with too low energy : %8d\n",
                        stopCount[6]);
                System.out.printf("              TOTAL RELEVANT STOPS : %8d\n",
                        stopCount[2] + stopCount[3] + stopCount[4] + stopCount[5]
                        + stopCount[6]);
                System.out.printf("\n");
            }

            // Fit residual plots
            for (int ri=0; ri<rn; ++ri) {
                double avgMean = 0;
                for (int li=1; li<=ln; ++li) {
                    avgMean += Math.abs(Data.fitRes(dgFMT[zi*rn + ri]
                            .getH1F("hi_cluster_res_l"+li),
                            dgFMT[zi*rn + ri].getF1D("f1_res_l"+li),
                            gssRanArr[ri]));
                }
                avgMeanArr[zi] = avgMean;
            }
        }

        // === DRAW PLOTS ======================================================
        for (int zi=0; zi<zn; ++zi) {
            System.out.printf("z shift : %4.1f, avg mean : %9.6f\n",
                    zShArr[zi], avgMeanArr[zi]);
        }

        EmbeddedCanvasTabbed fmtCanvas = Data.drawPlots(dgFMT, cn, titleArr, pltLArr);
        JFrame frame = new JFrame("FMT");
        frame.setSize(1600, 1000);
        frame.add(fmtCanvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);

        return 0;
    }
}
