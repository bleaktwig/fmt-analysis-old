/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.analysis;

import javax.swing.JFrame;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;

public class Data {
    /**
     * Fit residuals plot using two gaussians.
     * @param hires : 1-dimensional residuals distribution.
     * @param f1res : functions to use for fitting.
     * @param r     : radius for the gaussian fit (-r,r).
     * @return status int.
     */
    public static int fitRes(H1F hires, F1D f1res, double r) {
        // Fit the data
        double mean  = hires.getDataX(hires.getMaximumBin());
        double amp   = hires.getBinContent(hires.getMaximumBin());
        f1res.setParameter(0, amp);
        f1res.setParameter(1, mean);
        f1res.setParameter(2, 0.5);
        f1res.setRange(-r,r);
        DataFitter.fit(f1res, hires, "Q"); // No options uses error for sigma
	hires.setFunction(null);

        return 0;
    }

    /**
     * Get the mean of a distribution.
     * @param hires : 1-dimensional residuals distribution.
     * @return mean of hires.
     */
    public static double getMean(H1F hires) {
        double mean  = 0.0;
        double summ  = 0.0;
        double count = 0;
        for (int i = 0; i < hires.getAxis().getNBins(); ++i){
            summ  += hires.getAxis().getBinCenter(i) * hires.getBinContent(i);
            count += hires.getBinContent(i);
        }

        if (count!=0){
            mean = summ/count;
        }

        return mean;
    }

    /**
     * Get the sigma of a distribution.
     * @param hires : 1-dimensional residuals distribution.
     * @param mean  : mean of the distribution.
     * @return sigma of hires.
     */
    public static double getSigma(H1F hires, double mean) {
        double sigma = 0.0;
        double summ  = 0.0;
        double count = 0;
        for (int i=0; i < hires.getAxis().getNBins(); ++i) {
            double val = hires.getAxis().getBinCenter(i)*hires.getBinContent(i) - mean;
            summ  += val*val;
            count += hires.getBinContent(i);
        }
        if (count != 0) sigma = Math.sqrt(summ/(count-1));

        return sigma;
    }

    /**
     * Create a 1D data group of a specified variable.
     * @param var : index defining which variable is to be plotted.
     *                * 0 : cluster Tmin.
     *                * 1 : cluster total energy.
     * @param ln  : number of FMT layers.
     * @param r   : range of the variable (0,r).
     * @return the created data group.
     */
    public static DataGroup[] create1DDataGroup(int var, int ln, int r) {
        DataGroup[] dgFMT = new DataGroup[1];

        dgFMT[0] = new DataGroup(3,1);
        for (int li=1; li<=ln; ++li) {
            H1F hi_cluster_var =
                    new H1F("hi_cluster_var"+li, "", 200, 0, r);
            if (var == 0) hi_cluster_var.setTitleX("Tmin (ns)");
            if (var == 1) hi_cluster_var.setTitleX("energy");
            hi_cluster_var.setTitleY("cluster count");
            hi_cluster_var.setFillColor(4);

            dgFMT[0].addDataSet(hi_cluster_var, li-1);
        }

        return dgFMT;
    }

    /**
     * Create a 2D data group of a specified variable combination.
     * @param var : index defining which pair of variables is to be plotted:
     *                * 0 : energy vs cluster size.
     * @param ln  : number of FMT layers.
     * @param r   : unused.
     * @return the created data group.
     */
    public static DataGroup[] create2DDataGroup(int var, int ln, int r) {
        DataGroup[] dgFMT = new DataGroup[1];

        dgFMT[0] = new DataGroup(3,1);
        for (int li=1; li<=ln; ++li) {
            H2F hi_cluster_var = new H2F(
                    "hi_cluster_var"+li, "",
                    40, 0, 40, 100, 0, 5000
            );
            hi_cluster_var.setTitleX("cluster size");
            hi_cluster_var.setTitleY("energy / cluster size");
            dgFMT[0].addDataSet(hi_cluster_var, li-1);
        }

        return dgFMT;
    }

    /**
     * Render plots for a set of data groups.
     * @param dgFMT : Array of data groups.
     * @param title : Title to be given to the canvas.
     * @return status int.
     */
    public static int drawPlots(DataGroup[] dgFMT, String title) {
        EmbeddedCanvasTabbed fmtCanvas = new EmbeddedCanvasTabbed(title);

        fmtCanvas.getCanvas(title).draw(dgFMT[0]);
        fmtCanvas.getCanvas(title).setGridX(false);
        fmtCanvas.getCanvas(title).setGridY(false);
        fmtCanvas.getCanvas(title).setAxisFontSize(18);
        fmtCanvas.getCanvas(title).setAxisTitleSize(24);

        JFrame frame = new JFrame("FMT");
        frame.setSize(1600, 1000);
        frame.add(fmtCanvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);

        return 0;
    }

    /**
     * Create a set of data groups for residuals analysis.
     * @param ln : Number of FMT layers.
     * @param zn : Number of z shifts to try.
     * @param r  : Plot range.
     * @return array of data groups.
     */
    public static DataGroup[] createResDataGroups(int type, int ln, int zn, int r) {
        DataGroup[] dgFMT = new DataGroup[zn];

        for (int zi=0; zi<zn; ++zi) {
            dgFMT[zi] = new DataGroup(3,2);
            for (int li=1; li<=ln; ++li) {
                H1F hi_cluster_res =
                        new H1F("hi_cluster_res_l"+li, "", 200, -r, r);
                hi_cluster_res.setTitleX("Residual (cm) - Layer "+li);
                hi_cluster_res.setFillColor(4);
                dgFMT[zi].addDataSet(hi_cluster_res, li-1);

                F1D f1_res = new F1D(
                        "f1_res_l"+li,
                        "[amp]*gaus(x,[mean],[sigma])+[p0]+[p1]*x+[p2]*x*x",
                        -r,r
                );

                f1_res.setParameter(0, 0);
                f1_res.setParameter(1, 0);
                f1_res.setParameter(2, 1.0);
                f1_res.setLineWidth(2);
                f1_res.setLineColor(2);
                f1_res.setOptStat("1111");
                dgFMT[zi].addDataSet(f1_res, li-1);

                if (type == 0 || type == 1) {
                    H2F hi_cluster_res_strip = new H2F(
                        "hi_cluster_res_strip_l"+li,
                        200, -r, r, 100, 0, 1024
                    );
                    hi_cluster_res_strip.setTitleX("Residual (cm) - Layer "+li);
                    hi_cluster_res_strip.setTitleY("Strip - Layer "+li);
                    dgFMT[zi].addDataSet(hi_cluster_res_strip, li-1+ln);
                }
                if (type == 2) {
                    H2F hi_cluster_res_theta = new H2F(
                        "hi_cluster_res_theta_l"+li,
                        200, -r, r, 100, 0, 1
                    );
                    hi_cluster_res_theta.setTitleX("Residual - Layer "+li);
                    hi_cluster_res_theta.setTitleY("Theta - Layer "+li);
                    dgFMT[zi].addDataSet(hi_cluster_res_theta, li-1+ln);
                }
            }
        }

        return dgFMT;
    }

    /**
     * Render plots for a set of data groups, related to residuals analysis.
     * @param dgFMT    : Array containing the data groups.
     * @param cn       : Number of canvases to be drawn.
     * @param titleArr : Array of titles for the plots.
     * @param pltLArr  : Array of boolean defining which lines should be drawn
     *                   over the plots. Explained in detail in the
     *                   ResolutionAnalysis class.
     * @return status int.
     */
    public static int drawResPlots(DataGroup[] dgFMT, int cn,
            String[] titleArr, boolean[] pltLArr) {

        EmbeddedCanvasTabbed fmtCanvas = new EmbeddedCanvasTabbed(titleArr);
        for (int ci=0; ci<cn; ++ci) {
            fmtCanvas.getCanvas(titleArr[ci]).draw(dgFMT[ci]);
            fmtCanvas.getCanvas(titleArr[ci]).setGridX(false);
            fmtCanvas.getCanvas(titleArr[ci]).setGridY(false);
            fmtCanvas.getCanvas(titleArr[ci]).setAxisFontSize(18);
            fmtCanvas.getCanvas(titleArr[ci]).setAxisTitleSize(24);

            // Top plots
            for (int pi=0; pi<3; ++pi) {
                if (pltLArr[0]) {
                    DataLine vline =
                            new DataLine(0,0, 0,Double.POSITIVE_INFINITY);
                    vline.setLineColor(2);
                    vline.setLineWidth(2);
                    fmtCanvas.getCanvas(titleArr[ci]).cd(pi).draw(vline);
                }
            }

            // Bottom plots
            for (int pi=3; pi<6; ++pi) {
                if (pltLArr[1]) {
                    DataLine vline =
                            new DataLine(0,0, 0,Double.POSITIVE_INFINITY);
                    vline.setLineColor(0);
                    vline.setLineWidth(2);
                    fmtCanvas.getCanvas(titleArr[ci]).cd(pi).draw(vline);
                }
                if (pltLArr[2]) {
                    for (int j = 0; j < 16; ++j) {
                        DataLine hline =
                                new DataLine(-30,j*64 - 1, 30,j*64 - 1);
                        hline.setLineColor(0);
                        hline.setLineWidth(2);
                        fmtCanvas.getCanvas(titleArr[ci]).cd(pi).draw(hline);
                    }
                }

                if (pltLArr[3]) {
                    int[] seps = new int[]{320, 512, 832};
                    for (int sep : seps) {
                        DataLine hline = new DataLine(-30,sep, 30,sep);
                        hline.setLineColor(0);
                        hline.setLineWidth(2);
                        fmtCanvas.getCanvas(titleArr[ci]).cd(pi).draw(hline);
                    }
                }
            }
        }

        JFrame frame = new JFrame("FMT");
        frame.setSize(1600, 1000);
        frame.add(fmtCanvas);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);

        return 0;
    }
}
