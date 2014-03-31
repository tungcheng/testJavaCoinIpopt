/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package coinipopttest;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 *
 * @author TungNT
 */
public class CoinIpoptTest {

    /**
     * @param args the command line arguments
     */
    public static void main(String []args){
        
//        // Create the problem
//        TestBQP hs071 = new TestBQP();
//
//        // Get the default initial guess
//        double x[] = hs071.getInitialGuess();
//
//        // solve the problem
//        hs071.solve(x);
//        
//        for(int i=0; i<x.length; i++) {
//            System.out.println(x[i]);
//        }
        
        // Create the problem
        TestDCA test = new TestDCA();
//        TestBQP test = new TestBQP();

        // Get the default initial guess
        double x[] = test.getInitialGuess();

        // solve the problem
        test.solve(x);
        
        for(int i=0; i<x.length; i++) {
            System.out.println(x[i]);
        }
        
        
//        double[][] Q = {
//            { 0, 0, 0.5, 0.5},
//            { 0, 0, -0.5, 0.5},
//            { 0.5, -0.5, 0, 0.5},
//            { 0.5, 0.5, 0.5, 0}
//        };
//        double[][] x = {
//            { 1 },
//            { 1 },
//            { 1 },
//            { 1 }
//        };
//        Matrix mQ = new Matrix(Q);
//        mQ.print(4, 1);
//        Matrix mX = new Matrix(x);
//        mX.print(4, 1);
//        Matrix mt = getInverse(mX).times(mQ);
//        Matrix mt2 = mt.times(mX);
//        double value = mt2.get(0, 0);
//        System.out.println("value: " + value);
    }
    
    private static double getSmallestEigenvalue(Matrix Q) {
        EigenvalueDecomposition ed = Q.eig();
        Matrix x = ed.getD();
        double min = x.get(0, 0);
        for(int i=0; i < Q.getRowDimension(); i++) {
            if(min > x.get(i, i)) {
                min = x.get(i, i);
            }
        }
        return min;
    }
    
    private static Matrix getInverse(Matrix mX) {
        Matrix mXt = new Matrix( mX.getColumnDimension(), mX.getRowDimension());
        for(int i=0; i<mXt.getRowDimension(); i++) {
            for(int j=0; j<mXt.getColumnDimension(); j++) {
                mXt.set(i, j, mX.get(j, i));
            }
        }
        return mXt;
    }
}
