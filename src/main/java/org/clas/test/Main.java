/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.test;

import org.clas.analysis.TrkSwim;
import org.clas.analysis.ResolutionAnalysis;

public class Main {
        private static int usage() {
        System.out.printf("Usage: program infile\n");
        System.out.printf("  * infile: String with the hipo input file.\n");

        return 0;
    }

    public static void main(String[] args) {
        // Process input
        if (args.length != 1) {
            usage();
            System.exit(1);
        }

        String infile = args[0];

        // Setup
        double[] zShArr   = new double[]{-5, -4, -3};
        int pltRan        = 10;
        int gssRan        = 8;
        boolean[] pltLArr = new boolean[]{true, true, false, false};
        double[] swmSetup = new double[]{-0.75, -1.0, -3.0};
        boolean dbgInfo   = true;
        boolean testRun   = true;

        TrkSwim trkSwim = new TrkSwim(swmSetup);
        ResolutionAnalysis resAnls =
                new ResolutionAnalysis(infile, pltLArr, dbgInfo, testRun);

        // Run
        resAnls.zShiftAnalysis(zShArr, pltRan, gssRan, trkSwim);

        return;
    }
}
