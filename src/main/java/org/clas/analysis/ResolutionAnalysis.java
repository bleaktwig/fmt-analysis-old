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
     * Class constructor.
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

    /**
     * Print debugging information.
     * @param tsc : Cut tracks counter.
     * @param csc : Cut clusters counter.
     * @param en  : Total number of events read.
     * @return status int.
     */
    private int printDebugInfo(int[] tsc, int[] csc) {
        int tscsum = tsc[1] + tsc[2] + tsc[3] + tsc[4] + tsc[5];
        int cscsum = csc[1] + csc[2] + csc[3] + csc[4];
        System.out.printf("\n");
        System.out.printf("  tracks with wrong z coordinate │ %8d (%5.2f%%)   │\n", tsc[1], 100*((double)tsc[1])/tsc[0]);
        System.out.printf("     tracks not in the FMT layer │ %8d (%5.2f%%)   │\n", tsc[2], 100*((double)tsc[2])/tsc[0]);
        System.out.printf("     tracks in layer's bad areas │ %8d (%5.2f%%)   │\n", tsc[3], 100*((double)tsc[3])/tsc[0]);
        System.out.printf("     tracks with theta too large │ %8d (%5.2f%%)   │\n", tsc[4], 100*((double)tsc[4])/tsc[0]);
        System.out.printf(" tracks' with inconsistent Tmins │ %8d (%5.2f%%)   │\n", tsc[5], 100*((double)tsc[5])/tsc[0]);
        System.out.printf("            TOTAL TRACKS DROPPED │ %8d / %8d │\n", tscsum, tsc[0]);
        System.out.printf("                               %% │ %5.2f%%              │\n", 100*((double)tscsum)/tsc[0]);
        System.out.printf("─────────────────────────────────┼─────────────────────┤\n");
        System.out.printf("clusters with wrong strip number │ %8d (%5.2f%%)   │\n", csc[1], 100*((double)csc[1])/csc[0]);
        System.out.printf(" clusters with inappropiate Tmin │ %8d (%5.2f%%)   │\n", csc[2], 100*((double)csc[2])/csc[0]);
        System.out.printf("  small clusters with low energy │ %8d (%5.2f%%)   │\n", csc[3], 100*((double)csc[3])/csc[0]);
        System.out.printf("                clusters too big │ %8d (%5.2f%%)   │\n", csc[4], 100*((double)csc[4])/csc[0]);
        System.out.printf("          TOTAL CLUSTERS DROPPED │ %8d / %8d │\n", cscsum, csc[0]);
        System.out.printf("                               %% │ %5.2f%%              │\n", 100*((double)cscsum)/csc[0]);
        System.out.printf("─────────────────────────────────┴─────────────────────┘\n");

        return 0;
    }

    /**
     * Get data bank.
     * @param event : hipo event from which the bank is to be taken.
     * @param name  : data bank name.
     * @return reference to the data bank.
     */
    private DataBank getBank(DataEvent event, String name) {
        return (event.hasBank(name)) ? event.getBank(name) : null;
    }

    /**
     * Check if a track needs to be cut by z, x, y, or its theta angle. Tracks
     * that need to be cut due to the delta Tmin between clusters are processed
     * by another method due to its complexity.
     * @param z     : track's z coordinate at FMT layer.
     * @param x     : track's x coordinate at FMT layer.
     * @param y     : track's y coordinate at FMT layer
     * @param zRef  : FMT layer's z coordinate.
     * @param costh : cosine of track's theta angle (pz/p).
     * @param sc    : track error counter, used for debugging.
     * @return true if the track is to be cut, false otherwise.
     */
    private boolean checkTrackCuts(double z, double x, double y, double zRef,
                        double costh, int[] sc) {

        if (Math.abs(z-zRef)>0.05) {
            sc[2]++;
            return true;
        }
        if (25.0 > x*x + y*y || x*x + y*y > 225.0) {
            sc[3]++;
            return true;
        }
        if (costh>0.4) {
            sc[4]++;
            return true;
        }

        return false;
    }

    /**
     * Check if a cluster needs to be cut by its seed strip, size, energy, or
     * Tmin.
     * @param strip : seed strip of the cluster.
     * @param size  : cluster's size.
     * @param E     : cluster's total energy.
     * @param Tmin  : cluster's Tmin.
     * @param sc    : cluster error counter, used for debugging.
     * @return true if the cluster is to be cut, false otherwise.
     */
    private boolean checkClusterCuts(int strip, int size, double E, double Tmin,
            int[] sc) {

        if (strip<0 || strip>1023) {
            sc[1]++;
            return true;
        }
        if (Tmin < 50 || Tmin > 500) {
            sc[2]++;
            return true;
        }
        if (size == 1 && E < 100) {
            sc[3]++;
            return true;
        }
        if (size >= 5) {
            sc[4]++;
            return true;
        }

        return false;
    }

    /**
     * Generic function for running analysis, called by all others.
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
     * @return status int.
     */
    private int runAnalysis(int type, int opt, TrkSwim swim, DataGroup[] dgFMT,
            double zSh, int g) {
        // Sanitize input
        if (type < 0 || type > 2) return 1;

        int[] trackStopCount   = new int[]{0, 0, 0, 0, 0, 0};
        int[] clusterStopCount = new int[]{0, 0, 0, 0, 0};
        int ei = 0; // Event number.
        HipoDataSource reader = new HipoDataSource();
        reader.open(infile);
        System.out.printf("\nRunning analysis...\n");

        // Loop through events.
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

            if (type == 1 || type == 2) {
                track = getBank(event, "REC::Track");
                if (track==null) continue;
            }

            // Loop through trajectory points.
            for (int tri=0; tri<traj.rows(); tri++) {
                // Load trajectory variables.
                int detector = traj.getByte("detector", tri);
                int li = traj.getByte("layer", tri);
                int pi = traj.getShort("pindex", tri);
                int si = -1; // DC sector.
                double costh = -1; // track theta.

                // Use only FMT layers 1, 2, and 3.
                if (detector!=DetectorType.FMT.getDetectorId() || li<1 || li>ln)
                    continue;

                trackStopCount[0]++;

                // Get DC sector of the track.
                if (type == 1 || type == 2) {
                    for (int row = 0; row<track.rows(); ++row) {
                        if (track.getShort("pindex", row) == pi) {
                            si = track.getByte("sector", row)-1;
                        }
                    }
                }

                // Get FMT layer's z coordinate and phi angle.
                double zRef = fmtZ[li-1]/10 + zSh;
                double phiRef = fmtAngle[li-1];

                // Get particle data.
                double x  = (double) particle.getFloat("vx", pi);
                double y  = (double) particle.getFloat("vy", pi);
                double z  = (double) particle.getFloat("vz", pi);
                double px = (double) particle.getFloat("px", pi);
                double py = (double) particle.getFloat("py", pi);
                double pz = (double) particle.getFloat("pz", pi);
                int q     = (int)    particle.getByte("charge", pi);

                // Check if the track is too downstream before swimming.
                if (z > zRef) {
                    trackStopCount[1]++;
                    continue;
                }

                double[] V = swim.swimToPlane(x,y,z,px,py,pz,q,zRef);

                x  = V[0];
                y  = V[1];
                z  = V[2];
                px = V[3];
                py = V[4];
                pz = V[5];

                // Get the track's theta angle.
                costh = Math.acos(pz/Math.sqrt(px*px+py*py+pz*pz));

                // Apply track fiducial cuts.
                if (checkTrackCuts(z, x, y, zRef, costh, trackStopCount)) continue;

                // Rotate (x,y) to local coordinates.
                double xLoc = x * Math.cos(Math.toRadians(phiRef))
                        + y * Math.sin(Math.toRadians(phiRef));
                double yLoc = y * Math.cos(Math.toRadians(phiRef))
                        - x * Math.sin(Math.toRadians(phiRef));

                // Loop over the clusters and calculate residuals for every
                // track-cluster combination.
                for (int cri=0; cri<clusters.rows(); cri++) {
                    // Check that the cluster and trajectory layers match.
                    if (li!=clusters.getByte("layer", cri)) continue;
                    clusterStopCount[0]++;
                    int strip     = clusters.getInt("seedStrip", cri);
                    int size      = clusters.getShort("size", cri);
                    double energy = clusters.getFloat("ETot", cri);
                    double tmin   = clusters.getFloat("Tmin", cri);
                    double yclus  = clusters.getFloat("centroid", cri);

                    // Apply cluster fiducial cuts.
                    if (checkClusterCuts(strip, size, energy, tmin, clusterStopCount))
                        continue;

                    // Update plots depending on type of analysis.
                    int plti = -1;
                    if (type == 0) plti = opt;
                    if (type == 1 || type == 2) plti = si;

                    dgFMT[plti].getH1F("hi_cluster_res_l"+li)
                            .fill(yLoc - yclus);
                    if (type == 0 || type == 1)
                        dgFMT[plti].getH2F("hi_cluster_res_strip_l"+li)
                                .fill(yLoc - yclus, strip);
                    if (type == 2)
                        dgFMT[plti].getH2F("hi_cluster_res_theta_l"+li)
                                .fill(yLoc - yclus, costh);
                }
            }
        }
        System.out.format("Analyzed %8d events... Done!\n", ei);
        reader.close();

        if (debugInfo) printDebugInfo(trackStopCount, clusterStopCount);

        // Fit residual plots
        if (type == 0) {
            for (int li=1; li<=ln; ++li) {
                Data.fitRes(dgFMT[opt].getH1F("hi_cluster_res_l"+li),
                        dgFMT[opt].getF1D("f1_res_l"+li), g);
            }
        }
        else if (type == 1 || type == 2) {
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
     * Run z shift analysis.
     * @param zShArr  : Array of global z shifts to try.
     * @param r       : Plot range.
     * @param g       : Gaussian range.
     * @param swim    : TrkSwim class instance.
     * @return status int.
     */
    public int zShiftAnalysis(double[] zShArr, int r, int g, TrkSwim swim) {
        int type = 0;
        // === SETUP ===========================================================
        // Set canvases' stuff.
        int zn = zShArr.length;
        String[] titleArr = new String[zn];
        for (int zi=0; zi < zn; ++zi) titleArr[zi] = "z shift : "+zShArr[zi];

        DataGroup[] dgFMT = Data.createResDataGroups(type, ln, zn, r);

        double[] absMeanArr    = new double[zn];
        double[] absMeanErrArr = new double[zn];
        double[] sigmaArr     = new double[zn];
        double[] sigmaErrArr  = new double[zn];

        // === RUN =============================================================
        for (int zi=0; zi<zn; ++zi) {
            runAnalysis(type, zi, swim, dgFMT, zShArr[zi], g);

            // Get absolute mean and standard deviation from fit.
            for (int li=1; li<=ln; ++li) {
                F1D func = dgFMT[zi].getF1D("f1_res_l"+li);
                absMeanArr[zi]    += Math.abs(func.getParameter(1));
                absMeanErrArr[zi] += func.parameter(1).error();
                sigmaArr[zi]      += func.getParameter(2);
                sigmaErrArr[zi]   += func.parameter(2).error();
            }

            absMeanArr[zi]    /= ln;
            absMeanErrArr[zi] /= ln;
            sigmaArr[zi]     /= ln;
            sigmaErrArr[zi]  /= ln;
        }

        // === PRINT ALIGNMENT DATA AND DRAW PLOTS =============================
        for (int zi=0; zi<zn; ++zi) {
            System.out.printf("z shift : %5.2f\n", zShArr[zi]);
            System.out.printf("  * |mean| : %9.6f +- %9.6f\n",
                    absMeanArr[zi], absMeanErrArr[zi]);
            System.out.printf("  * sigma  : %9.6f +- %9.6f\n",
                    sigmaArr[zi], sigmaErrArr[zi]);
        }

        Data.drawResPlots(dgFMT, zn, titleArr, pltLArr);

        return 0;
    }

    /**
     * Run DC sector strip analysis.
     * @param zSh  : z shift.
     * @param r    : Plot range.
     * @param g    : Gaussian range.
     * @param swim : TrkSwim class instance.
     * @return status int.
     */
    public int dcSectorStripAnalysis(double zSh, int r, int g, TrkSwim swim) {
        int type = 1;
        // Set canvases' stuff.
        int sn = 6; // Number of DC sectors.
        String[] titleArr = new String[sn];
        for (int si=0; si < sn; ++si) titleArr[si] = "DC sector "+(si+1);

        // Run.
        DataGroup[] dgFMT = Data.createResDataGroups(type, ln, sn, r);
        runAnalysis(type, sn, swim, dgFMT, zSh, g);
        Data.drawResPlots(dgFMT, sn, titleArr, pltLArr);

        return 0;
    }

    /**
     * Run DC sector theta analysis.
     * @param zSh  : z shift.
     * @param r    : Plot range.
     * @param g    : Gaussian range.
     * @param swim : TrkSwim class instance.
     * @return status int.
     */
    public int dcSectorThetaAnalysis(double zSh, int r, int g, TrkSwim swim) {
        int type = 2;
        // Set canvases' stuff.
        int sn = 6; // Number of DC sectors.
        String[] titleArr = new String[sn];
        for (int si=0; si < sn; ++si) titleArr[si] = "DC sector "+(si+1);

        // Run.
        DataGroup[] dgFMT = Data.createResDataGroups(type, ln, sn, r);
        runAnalysis(type, sn, swim, dgFMT, zSh, g);
        Data.drawResPlots(dgFMT, sn, titleArr, pltLArr);

        return 0;
    }

    /**
     * Draw a 1D plot by counting a pre-defined variable.
     * @param var : Variable to be counted:
     *              * 0 : clusters' Tmin.
     *              * 1 : clusters' energy.
     * @param r   : Range for the plot (min = 0, max = r)
     * @return status int.
     */
    public int plot1DCount(int var, int r) {
        // Sanitize input
        if (var < 0 || var > 1) return 0;

        String title = null;
        if (var == 0) title = "Tmin count";
        if (var == 1) title = "energy count";

        DataGroup[] dgFMT = Data.create1DDataGroup(var, ln, r);

        // Run
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
            if (clusters==null || traj==null) continue;

            for (int ri=0; ri<clusters.rows(); ++ri) {
                int li        = clusters.getByte("layer", ri);
                double energy = clusters.getFloat("ETot", ri);
                double tmin   = clusters.getFloat("Tmin", ri);

                if (var==0) dgFMT[0].getH1F("hi_cluster_var"+li).fill(tmin);
                if (var==1) dgFMT[0].getH1F("hi_cluster_var"+li).fill(energy);
            }
        }
        System.out.format("Analyzed %8d events... Done!\n", ei);
        reader.close();

        Data.drawPlots(dgFMT, title);

        return 0;
    }

    /**
     * Draw a 2D plot by counting a pre-defined variable against another.
     * @param var : Variable pair to be counted:
     *              * 0 : energy / cluster size vs cluster size.
     * @param r   : Currently unused, kept for consistency.
     * @return status int.
     */
    public int plot2DCount(int var, int r) {
        // Sanitize input
        if (var != 0) return 0;

        String title = null;
        if (var == 0) title = "energy / cluster size count";

        DataGroup[] dgFMT = Data.create2DDataGroup(var, ln, r);

        // Run
        int ei = 0; // Event number.
        HipoDataSource reader = new HipoDataSource();
        reader.open(infile);
        System.out.printf("\nRunning analysis...\n");

        // === LOOP THROUGH EVENTS =============================================
        while (reader.hasEvent()) {
            if (ei == 10000 && testRun) break;
            if (ei%50000==0) System.out.format("Analyzed %8d events...\n", ei);
            DataEvent event = reader.getNextEvent();
            ei++;

            // Get relevant data banks.
            DataBank clusters = getBank(event, "FMTRec::Clusters");
            DataBank traj     = getBank(event, "REC::Traj");
            if (clusters==null || traj==null) continue;

            for (int ri=0; ri<clusters.rows(); ++ri) {
                int li        = clusters.getByte("layer", ri);
                double energy = clusters.getFloat("ETot", ri);
                double tmin   = clusters.getFloat("Tmin", ri);
                int size      = clusters.getShort("size", ri);

                // if (applyCuts(energy, tmin)) continue;

                if (var==0) dgFMT[0].getH2F("hi_cluster_var"+li)
                                    .fill(size, energy/size);
            }
        }
        System.out.format("Analyzed %8d events... Done!\n", ei);
        reader.close();

        Data.drawPlots(dgFMT, title);

        return 0;
    }
}
