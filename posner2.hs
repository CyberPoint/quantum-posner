{-

author: James Ulrich julrich@cyberpointllc.com


Copyright (c) 2016 CyberPoint International LLC


Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
-}

import Quipper
import Data.List.Split
import QuipperLib.Synthesis
import Quantum.Synthesis.Ring
import Quantum.Synthesis.Matrix
import Quantum.Synthesis.MultiQubitSynthesis
import Debug.Trace
import QuipperLib.Simulation
import System.Random
import Data.Map (Map)
import qualified Data.Map as Map
import Quipper.QData

{-

This program will take as input a 1-1 map of key-value pairs, where the key is 
an n-digit binary number (an "input"), and the value is an n-digit binary 
number (the corresponding "output"). It will produce as output a quantum 
circuit, which with high probability, maps a given input to the corresponding 
output.  

We assume program is invoked as 'posner2 < func.txt' where func.txt is a file 
of the  form (here using example n=2):
        
    00:00   
    01:01
    10:10
    11:11
       01
   <mode>
                   
the first 2^n lines define a bijective function from bitstrings of length n to 
bitstrings of length n; the digits to the left of the ':' are the inputs. the 
lines must appear in increasing order of  outputs, viewed as binary strings 
(e.g. 00, 01, 10, 11, for n=4). the next to last line 
is a specific input to be evaluated by simulating the circuit (used for mode "S"
only). the last line is a run mode: 
    
    A=show amplitudes via sim_amps (supported for n <= 2 only), for input |0>.
    P=print_generic.
    S=simulate via run_generic.


-- major algorithmic steps (in concept, if not actual code) are:

1. ingest a string giving the map of (input,output) pairs and 
   store in a suitable data structure in memory.
   
2. define an n-quibit circuit that takes the |0> state to 
   one that is a superposition S of all possible n-qubit basis states, 
   but for which the state |0> has the largest probability amplitude. We do
   by first applying a hadamard gate, and then applying the two level "bias" 
   gate  [[cos(pi/6), sin(pi/6)], [sin(pi/6), -cos(pi/6)]] successively to 
   states (|0>,|1>), (|0>,|2>), ...., (|0>,|2^(n-1)>). 
   
   so if you consider what happens to amplitude of |0> under this 
   transformation, starting with amplitude of 1:
   
   first it goes to a_1 = (1/sqrt(2)^n)[cos(pi/6) + sin(pi/6)].
   then to a_2 = a_1 cos(pi/6) + (1/sqrt(2)^n)sin(pi/6).
   then to a_3 = a_2 cos(pi/6) + (1/sqrt(2)^n)sin(pi/6)....
   
   so it ends up as a_(n-1) cos(pi/6) + (1/sqrt(2)^n) sin(pi/6)
   where a_i = a_{i-1} cos(pi/6) + (1/sqrt(2)^n) sin(pi/6). should be
   able to calculate that in python.  
     
   
3. add n ancilla qubits to the circuit. then for each output i given by the 
   map, add a controlled swap gate exchanging the weight of the |0> component 
   of S with that of the component corresponding to |i>, with controls given 
   by the input producing |i>. the controls are placed on the ancillas.
   
4. to evaluate the circuit for a given input, the ancillas are initialized to
   the input.
   
   see also http://arxiv.org/abs/1601.07137
    
-}

-- declare a map data type entry
data InputOutput = InputOutput {
   
   input :: [Bool],
   output :: [Bool]
} deriving (Show)

-- declare a function truth table (input->output map defining a function)
data TruthTable = TruthTable {
    entries :: [InputOutput]
} deriving (Show)
 

-- produce a quantum circuit from the function defined by the input 
processInput :: [String] -> ([Qubit] -> Circ [Qubit])
processInput table = do

    -- from input lines defining a function (as a text-based truth table of 
    -- inputs->outputs), get list of input,output string pairs
    let u = map (\x -> splitOn ":" x) table
    
    -- from list of input, output string pairs, get list of 
    -- input, output boolean list pairs, make internal truth table from that
    let v2 = map (\x -> (str_to_blist(x !! 0),str_to_blist(x !! 1)) ) u  
    let g2 = map (\x -> InputOutput (fst x) (snd x)) v2      
    let tt = TruthTable g2   

    -- make quantum circuit based on table 
    let circ = tt_circuit tt
    circ
    

-- convert a text string of 0s and 1s to a boolean array   
str_to_blist :: String -> [Bool]
str_to_blist lc = 
    let bl = map (\x -> case x of '1' -> True
                                  '0' -> False) lc
    in bl                             
                            

--  A basis change to obtain y-rotations from z-rotations. 
--  courtesy Neil J. Ross.                    
y_to_Z :: Qubit -> Circ Qubit
y_to_Z q = do
  q <- gate_S_inv q
  q <- gate_H q
  q <- gate_S_inv q
  return q


{-
the "bias" gate [ [cos(pi/6), sin(pi/6)], [sin(pi/6), -cos(pi.6)]] expressed as 
         
          U = e^{i \alpha) R_Z(\beta) R_Y(\gamma) R_Z(\delta) 
       
where \alpha = pi/2, \beta = 0, \delta = pi, \gamma = pi/3, which 
reduces to R_Y(\gamma)*Z; see Nielsen and Chuang, Quantum Computation and 
Quantum Information, pp.175-176. Decomposition in terms of Quipper gates 
courtesy Neil J. Ross.
-}
my_biasUt :: Timestep -> Qubit -> Circ ()
my_biasUt t q = do
  
    -- q <- expZt pi q
    q <- gate_Z q
    q <- with_computed (y_to_Z q) (expZt t)
    return ()


apply_controlled_X :: ControlList -> Qubit -> Circ()
apply_controlled_X controls q = do
    q <- gate_X q `controlled` controls 
    return ()
    
-- repeatedly apply our two-level bias gate  [[ sqrt(p), 0.5], [0.5, -sqrt(p)]] as 
-- a U_{0,i} gate for i, i-1, i-2, ..., 1
bias_gate_recurse  :: (Qubit-> Circ()) -> Int -> [Qubit] -> Circ()
bias_gate_recurse g 0 qbits =  return()
bias_gate_recurse g i qbits = do 
    twolevel 0 i qbits g
    bias_gate_recurse g (i-1) qbits  

-- construct a control list for a given set of function inputs, expressed 
-- as bools
control_recurse :: ControlList -> Int -> [Qubit] -> [Bool] -> ControlList
control_recurse controls idx qubits bools = 
    case (idx == 0) of 
        True -> controls
        False -> do 
            let new_controls = controls .&&. (qubits !! idx) .==. (bools !! idx)    
            control_recurse new_controls (idx-1) qubits bools

-- apply a controlled X gate as U_{0,idx} with controls on inputs corresponding
-- to output state |idx> for thruth table entries idx = 1....n                  
flip_gate_recurse ::  TruthTable -> Int -> [Qubit] -> [Qubit] -> Circ()
flip_gate_recurse tt idx ancillas qbits = do

     
    -- get the input booleans for the output given by idx
    let pos_neg = input ( (entries tt) !! idx ) 

    -- now build a control list based on the booleans. we recurse to build it 
    -- up wire by wire using the .&&. operator    
    let init_controls = (ancillas !! 0) .==. (pos_neg !! 0)    
    let controls = control_recurse init_controls (length(pos_neg) -1 ) 
                                    ancillas pos_neg
        
    let g = apply_controlled_X controls
    twolevel 0 idx qbits g

    case (idx == (length (entries tt) - 1)) of 
        True -> return()
        False -> flip_gate_recurse tt (idx + 1) ancillas qbits

       
-- construct a circuit that computes the output for the given input and truth 
-- table
tt_circuit :: TruthTable -> [Qubit] -> Circ [Qubit]
tt_circuit tt ancilla = do

    -- n = number of input bits
    let n = length (input ( (entries tt) !! 0) )
    
    -- create n circuit bits and set them all to 0s                                
    qbits <- qinit (replicate n False)
    
    -- now perform a hadamard on the qbits
    mapUnary hadamard qbits
   
    -- apply our 2-level bias gate to (|0>,|1>), (|0>,|2>), ..., (|0>,|2^(n-1)>)
    -- note:  we really want a y-rotation of pi/divisor but for some reason in quipper,
    -- we need to halve it (pi/2*divisor)
    --let divisor = if (even n) then 2*(2^((div n 2)+1)) else ( + 2*(2^((div n-1 2)+1)) 0.5*(2^((div n-1 2)+1)))
         
    let p = pi/(get_divisor n)
    let g = my_biasUt p  
    
    -- step 2
    bias_gate_recurse g (2^n -1) qbits
    
    -- step 3
    flip_gate_recurse tt 1 ancilla qbits
  
    qdiscard ancilla
    return qbits

get_divisor :: Int -> Double
get_divisor n = case (even n) of 
                    True -> 2*(2^((div n 2) + 1)) 
                    False ->  2*(2^((div (n-1) 2) + 1) +  (  0.5*(2^((div (n-1) 2) + 1))  ) ) 
                     
    
-- here follow some utility routines to help with debugging    
tester :: [Bool]  -> ([Qubit] -> Circ [Qubit]) -> Circ [Qubit]
tester ibools circ = do
    ibits <- qinit ibools 
    circ ibits
    
-- a sample input state for a 1-qbit circuit
myMap :: Map [Bool] (Cplx Double)
myMap = Map.fromList (map makePair [False, True])
    where makePair x =  case x of False -> ([False],1.0)
                                  True -> ([True], 0.0)

-- a sample input state for a 2-qbit circuit
myMap2 :: Map [Bool] (Cplx Double)
myMap2 = Map.fromList (map makePair
    [[False, False],[False,True], [True,False],[True,True]])
    where makePair x =  case x of [False,False] -> ([True, True],1.0)
                                  [False, True] -> ([False,True], 0.0)
                                  [True, False] -> ([True,False], 0.0)
                                  [True, True]  -> ([False,False],  0.0)
                                                                          
-- and here's main!                                        
main = do

    s <- getContents  
    let t = lines s 
    
    -- separate the lines defining the function to be evaluated from the 
    -- input on which it is to be evaluated, and the mode
    let func_def = take ((length t)- 2) t
    let input = t !!  ((length t) - 2)
    let mode = t !! ((length t) - 1)
    let ibools = str_to_blist input
    let ibits = map (\x -> cinit x) ibools
    
    -- get a quantum circuit that evaluates the function given an input 
    let qa = processInput func_def

    case mode of
        "A" -> do
            -- WE USE THIS CODE TO GET OUTPUT STATE AMPLITUDES FOR THE 
            -- FUNCTION DEFINED BY STDIN, AND THE INPUT HARDCODED IN myMap, 
            -- myMap2  RTNS 
            case ((length ibools) < 3) of 
                True -> do
                    let m = if (length ibools) == 1 then myMap else myMap2  
                    t <- randomIO
                    putStrLn("here is the input map: " ++  Map.showTree(m))      
                    let res2 = sim_amps (mkStdGen t) qa m     
                    let s3 = Map.showTree(res2)      
                    putStrLn("here is the output map: ")
                    putStrLn(s3)
                False -> do     
                    putStrLn("option A only supported for n=2")

        "P" -> do
            -- WE USE THIS CODE TO PRINT THE CIRCUIT, GIVEN THE FUNCTION
            -- DEFINED BY STDIN
            let circ = tester ibools qa 
            print_generic Preview circ
        
        "S" -> do   

            -- WE USE THIS CODE TO RUN MANY SIMULATION TRIALS OVER THE CIRCUIT,
            -- GIVEN THE FUNCTION AND INPUT DEFINED BY STDIN
            for 1 100 1 $ \i -> do  
                t <- randomIO
                let res = run_generic (mkStdGen t) (1.0::Double) qa ibools
                putStrLn(show(res))
            endfor

   