/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package coinipopttest;

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
//        HS071 hs071 = new HS071();
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
        TestBQP hs071 = new TestBQP();

        // Get the default initial guess
        double x[] = hs071.getInitialGuess();

        // solve the problem
        hs071.solve(x);
        
        for(int i=0; i<x.length; i++) {
            System.out.println(x[i]);
        }
        
        TestBQP hs071_2 = new TestBQP();
        // Get the default initial guess
        double x2[] = hs071_2.getInitialGuess();

        // solve the problem
        hs071_2.solve(x2);
        
        for(int i=0; i<x2.length; i++) {
            System.out.println(x2[i]);
        }
    }
    
}
