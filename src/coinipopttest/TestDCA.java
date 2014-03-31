/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package coinipopttest;

import Jama.Matrix;
import org.coinor.Ipopt;

/**
 *
 * @author TungNT
 */
public class TestDCA extends Ipopt {
    
    static final int T = 1000000;
    
    double[][] temp;
    
    double[][] Q = {
            { 27.8, -14.0, -18.0},
            { -14.0,  43.8, -22.0},
            { -18.0, -22.0,  37.8}
        };
    
    double[][] q = {
            { -26.9 },
            { -26.9 },
            { -26.9 }
        };
    
    Matrix mX, mQ, mq;
    
    // Problem sizes
    int n, m, nele_jac, nele_hess;
    
    public TestDCA() {
        n = 3;
        temp = new double[this.n][1];
        mX = new Matrix(temp);
        mQ = new Matrix(Q);
        mq = new Matrix(q);
        
        /* Number of nonzeros in the Jacobian of the constraints */
        nele_jac = 3;
        /* Number of nonzeros in the Hessian of the Lagrangian (lower or
         * upper triangual part only) */
        nele_hess = 6;

        /* set the number of variables and allocate space for the bounds */
        n=3;
        double x_L[] = new double[n];
        double x_U[] = new double[n];
        for(int i=0; i < x_L.length; i++){
                x_L[i] = 0.0;
                x_U[i] = 1.0;
        }

        /* set the number of constraints and allocate space for the bounds */
        m=1;
        double g_L[] = new double[m];
        double g_U[] = new double[m];
        /* set the values of the constraint bounds */
        g_U[0] = 2;
        g_L[0] = (-1)*2e19;

        /* Index style for the irow/jcol elements */
        int index_style = Ipopt.C_STYLE;

        /* create the IpoptProblem */
        create(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, index_style);
    }
    
    public double[] getInitialGuess(){
        /* allocate space for the initial point and set the values */
        double x[] = new double[3];
        x[0] = 0.5;
        x[1] = 0.5;
        x[2] = 0.5;

        return x;
    }

    protected boolean eval_f(int n, double[] x, boolean new_x, double[] obj_value) {
        assert n == this.n;
        
        for(int i=0; i<this.n; i++) {
            temp[i][0] = x[i];
        }

        //obj_value[0] = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
        Matrix t = mX.inverse().times(mQ);
        Matrix t2 = t.times(mX);
        Matrix t3 = mq.inverse().times(mX);
        
        obj_value[0] = 0.5 * t2.get(0, 0) + t3.get(0, 0);
        
        for(int i=0; i<this.n; i++) {
            obj_value[0] += T*Math.min(x[i], 1 - x[i]);
        }

        return true;
    }

    protected boolean eval_grad_f(int n, double[] x, boolean new_x, double[] grad_f) {
        assert n == this.n;

        for(int i=0; i<this.n; i++) {
            temp[i][0] = x[i];
        }
        
        Matrix t = mQ.inverse().plus(mQ);
        Matrix t2 = t.times(mX);
        Matrix t3 = t2.times(0.5);
        Matrix t4 = mq.plus(t3);
        
        for(int i=0; i<this.n; i++) {
            grad_f[i] = t4.get(i, 0);
            if(x[i] < 0.5) {
                grad_f[i] += T;
            }
            else {
                grad_f[i] += (-T);
            }
        }

        return true;
    }

    protected boolean eval_g(int n, double[] x, boolean new_x, int m, double[] g) {
            assert n == this.n;
            assert m == this.m;
            
            g[0] = 0;
            for(int i=0; i<n; i++) {
                g[0] += x[i];
            }

            return true;
    }

    protected boolean eval_jac_g(int n, double[] x, boolean new_x,
                    int m, int nele_jac, int[] iRow, int[] jCol, double[] values) {
            assert n == this.n;
            assert m == this.m;

            if (values == null) {
                    /* return the structure of the jacobian */

                    /* this particular jacobian is dense */
                    iRow[0] = 0;
                    jCol[0] = 0;
                    iRow[1] = 0;
                    jCol[1] = 1;
                    iRow[2] = 0;
                    jCol[2] = 2;
            }
            else {
                    /* return the values of the jacobian of the constraints */

                    values[0] = 1;
                    values[1] = 1;
                    values[2] = 1;
            }

            return true;
    }

    protected boolean eval_h(int n, double[] x, boolean new_x, double obj_factor, int m, double[] lambda, boolean new_lambda, int nele_hess, int[] iRow, int[] jCol, double[] values) {
            int idx = 0; /* nonzero element counter */
            int row = 0; /* row counter for loop */
            int col = 0; /* col counter for loop */
            if (values == null) {
                    /* return the structure. This is a symmetric matrix, fill the lower left
                     * triangle only. */

                    /* the hessian for this problem is actually dense */
                    idx=0;
                    for (row = 0; row < n; row++) {
                            for (col = 0; col <= row; col++) {
                                    iRow[idx] = row;
                                    jCol[idx] = col;
                                    idx++;
                            }
                    }

                    assert idx == nele_hess;
                    assert nele_hess == this.nele_hess;
            }
            else {
                Matrix t = mQ.inverse().plus(mQ);
                Matrix t2 = t.times(0.5);
                    /* return the values. This is a symmetric matrix, fill the lower left
                     * triangle only */

                    /* fill the objective portion */
                    values[0] = obj_factor * t2.get(0, 0);                 /* 0,0 */

                    values[1] = obj_factor * t2.get(1, 0);                 /* 1,0 */
                    values[2] = obj_factor * t2.get(1, 1);                 /* 1,1 */

                    values[3] = obj_factor * t2.get(2, 0);                 /* 2,0 */
                    values[4] = obj_factor * t2.get(2, 1);                 /* 2,1 */
                    values[5] = obj_factor * t2.get(2, 2);                 /* 2,2 */

            }
            return true;
    }    
}
