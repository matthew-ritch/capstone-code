/*
00002  *  Copyright 2006-2007 Columbia University.
00003  *
00004  *  This file is part of MEAPsoft.
00005  *
00006  *  MEAPsoft is free software; you can redistribute it and/or modify
00007  *  it under the terms of the GNU General Public License version 2 as
00008  *  published by the Free Software Foundation.
00009  *
00010  *  MEAPsoft is distributed in the hope that it will be useful, but
00011  *  WITHOUT ANY WARRANTY; without even the implied warranty of
00012  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
00013  *  General Public License for more details.
00014  *
00015  *  You should have received a copy of the GNU General Public License
00016  *  along with MEAPsoft; if not, write to the Free Software
00017  *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
00018  *  02110-1301 USA
00019  *
00020  *  See the file "COPYING" for the text of the license.
00021  */

package capstone.wheelchairmonitor;

public class FFT {

    int n, m;

    // Lookup tables. Only need to recompute when size of FFT changes.
    double[] cos;
    double[] sin;

    public FFT(int n) {
        this.n = n;
        this.m = (int) (Math.log(n) / Math.log(2));

        // Make sure n is a power of 2
        if (n != (1 << m))
            throw new RuntimeException("FFT length must be power of 2");

        // precompute tables
        cos = new double[n / 2];
        sin = new double[n / 2];

        for (int i = 0; i < n / 2; i++) {
            cos[i] = Math.cos(-2 * Math.PI * i / n);
            sin[i] = Math.sin(-2 * Math.PI * i / n);
        }

    }

    public void fft(double[] x, double[] y) {
        int i, j, k, n1, n2, a;
        double c, s, t1, t2;

        // Bit-reverse
        j = 0;
        n2 = n / 2;
        for (i = 1; i < n - 1; i++) {
            n1 = n2;
            while (j >= n1) {
                j = j - n1;
                n1 = n1 / 2;
            }
            j = j + n1;

            if (i < j) {
                t1 = x[i];
                x[i] = x[j];
                x[j] = t1;
                t1 = y[i];
                y[i] = y[j];
                y[j] = t1;
            }
        }

        // FFT
        n1 = 0;
        n2 = 1;

        for (i = 0; i < m; i++) {
            n1 = n2;
            n2 = n2 + n2;
            a = 0;

            for (j = 0; j < n1; j++) {
                c = cos[a];
                s = sin[a];
                a += 1 << (m - i - 1);

                for (k = j; k < n; k = k + n2) {
                    t1 = c * x[k + n1] - s * y[k + n1];
                    t2 = s * x[k + n1] + c * y[k + n1];
                    x[k + n1] = x[k] - t1;
                    y[k + n1] = y[k] - t2;
                    x[k] = x[k] + t1;
                    y[k] = y[k] + t2;
                }
            }
        }
    }
}