import java.util.*;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import java.util.Random;

public class LearningAutomata
{
    private static Random rnd = new Random();
    private static Random rndDouble = new Random();
    private static double alpha = .89;
    private static double beta = .19;
    
    private static short maxIndex(DoubleArrayList arrayList) {
        short index = 0;
        for(short i = 1; i < arrayList.size(); i++) {
            if(arrayList.get(i) > arrayList.get(index)) {
                index = i;
            }
        }
        return index;
    }
    
    private static double maxValue(DoubleArrayList arrayList) {
        return arrayList.get(maxIndex(arrayList));
    }
    
    private static double sum(DoubleArrayList arrayList) {
        double sum = 0;
        for(short i = 0; i < arrayList.size(); i++) {
            sum += arrayList.get(i);
        }
        return sum;
    }
    
    private static short sum(ShortArrayList arrayList) {
        short sum = 0;
        for(short i = 0; i < arrayList.size(); i++) {
            sum += arrayList.get(i);
        }
        return sum;
    }
    

    private static short rouletteWheel(DoubleArrayList probability) {
        short returnedIndex = -1;
        short numEvents = (short) probability.size();
        double maxProbability = maxValue(probability);
        double epsilon = 1e-6;
        short numberOfRandomNumbers = (short) (2 * numEvents);
        
        if((1 - maxProbability) < epsilon) {
            returnedIndex = maxIndex(probability);
        } else {
            ShortArrayList sequence = new ShortArrayList();
            for(short i = 0; i < numEvents; i++) {
                if(probability.get(i) > 0) {
                    short times = (short) (numberOfRandomNumbers * probability.get(i));
                    for(short j = 0; j < times; j++) {
                        sequence.add(i);
                    }
                }
            }
            Collections.shuffle(sequence);
            short randomNumber = (short) rnd.nextInt(sequence.size());
            returnedIndex = sequence.get(randomNumber);
            System.out.println(sequence);
            
        }
        return(returnedIndex);
    }

    
    private static short actionSelection(DoubleArrayList probability) {
        short returnedIndex = -1;
        short numEvents = (short) probability.size();
        double maxProbability = maxValue(probability);
        double epsilon = 1e-6;
        
        if((1 - maxProbability) < epsilon) {
            returnedIndex = maxIndex(probability);
        }
        else {
            //double[] tempProbabilities = new double[numEvenets];
            //for(int i = 0; i < numEvenets; i++) {
            //    tempProbabilities[i] = probability.get(i);
            //}
            DoubleArrayList tempProbabilities = new DoubleArrayList(probability);
            ShortArrayList tempIndices = new ShortArrayList();
            
            //short[] tempIndices = new short[numEvenets];
            for(short i = 0; i < numEvents; i++) {
                tempIndices.add(i);
            }
            
            short factor = 2;
            double seprator = (double) 1/factor;
            double randomNumber;
            
            while(numEvents > factor) {
                //double sumTempProbabilities = sum(tempProbabilities);
                //if(sumTempProbabilities != 1) {
                //    tempProbabilities.set(0, (double) (tempProbabilities.get(0) + 1 - sumTempProbabilities));
                //}
                
                
                randomNumber = rnd.nextDouble();
                //System.out.println("RND " + randomNumber);
                
                short leftMax = 0;
                short leftNum = 0;
                short rightMin = 0;
                short rightNum = 0;
                short minIndex = tempIndices.get(0);
                double sumProbabilities = 0;
                short k = 0;
                short j = 0;
                
                do 
                {
                    sumProbabilities += tempProbabilities.get(k);
                    //System.out.println(k + " " + sumProbabilities + " " + (sumProbabilities - seprator));
                    k++;
                } while((Math.abs(sumProbabilities - seprator) > epsilon) && (sumProbabilities < seprator));
                //System.out.println("SUM " + sumProbabilities);    
                
                rightMin = k;
                leftMax = (short) (rightMin - 1);
                leftNum = (short) (leftMax + 1);
                rightNum = (short) (numEvents - leftNum);
                if(Math.abs(sumProbabilities - seprator) > epsilon) {
                    rightMin--;
                    rightNum++;
                }	
                //System.out.println(leftMax + " " + leftNum + " " + rightMin + " " + rightNum);    
                if(randomNumber < seprator) {
                    
                    //System.out.println("Left ");
                    //short[] leftIndices = new short[leftNum];
                    DoubleArrayList leftProbabilities = new DoubleArrayList();
                    ShortArrayList leftIndices = new ShortArrayList();
                    
                    for(k = leftMax; k >= leftNum - leftMax - 1; k--) {
                        j = (short) (leftNum - k - 1);
                        leftProbabilities.add(tempProbabilities.get(j));
                        leftIndices.add((short) (j + minIndex));
                    }
                   // System.out.println("leftProbabilities" + leftProbabilities + " " + leftProbabilities.get(leftNum - 1) + " " + sum(leftProbabilities) + " " + seprator );   
                    if(Math.abs(sumProbabilities - seprator) > epsilon) {
                        leftProbabilities.set(leftNum - 1, leftProbabilities.get(leftNum - 1) - (sum(leftProbabilities) - seprator));
                    }
                    //System.out.println("leftProbabilities" + leftProbabilities);   
                    for(k = 0; k < leftNum; k++) {   
                        leftProbabilities.set(k, leftProbabilities.get(k) * factor);
                    }
                    //System.out.println("leftProbabilities" + leftProbabilities);   
                    //System.out.println("leftIndices" + leftIndices);   
                    
                    tempProbabilities = new DoubleArrayList(leftProbabilities);
                    tempIndices = new ShortArrayList(leftIndices);
                    numEvents = leftNum;                    
                }
                else {
                    /*
                    System.out.println("Right ");
                    System.out.println("tempProbabilities" + tempProbabilities); 
                    System.out.println("tempIndices" + tempIndices);
                    System.out.println("minIndex " + minIndex);
                    System.out.println("rightMin " + rightMin);
                    System.out.println("rightNum " + rightNum);
                    */
                    DoubleArrayList rightProbabilities = new DoubleArrayList();
                    ShortArrayList rightIndices = new ShortArrayList();
                    
                    for(k = rightMin; k < rightNum + rightMin; k++) {
                        j = (short) (k - rightMin);
                        rightProbabilities.add(tempProbabilities.get(k));
                        rightIndices.add((short) (k + minIndex));
                    }
                    //System.out.println("1rightProbabilities" + rightProbabilities);   
                    //System.out.println("1rightIndices" + rightIndices);   
                    if(Math.abs(sumProbabilities - seprator) > epsilon) {
                        rightProbabilities.set(0, rightProbabilities.get(0) - (sum(rightProbabilities) - seprator));
                    }
                    //System.out.println("2rightProbabilities" + rightProbabilities);   
                    for(k = 0; k < rightNum; k++) {   
                        rightProbabilities.set(k, rightProbabilities.get(k) * factor);                    
                    }
                    //System.out.println("3rightProbabilities" + rightProbabilities);   
                    tempProbabilities = new DoubleArrayList(rightProbabilities);
                    tempIndices = new ShortArrayList(rightIndices);
                    numEvents = rightNum;
                }
            }
            //System.out.println(tempProbabilities);
            //System.out.println(tempIndices);
            
            if(numEvents == 1) {
                returnedIndex = tempIndices.get(0); // Need to check this
            }
            else if(numEvents == factor)
            {
                randomNumber = rnd.nextDouble();
                if(randomNumber < tempProbabilities.get(0))
                    returnedIndex = tempIndices.get(0);
                else
                    returnedIndex = tempIndices.get(1);
            }
        }
        return returnedIndex;
    }
    
    private static void weightedProbabilityUpdate(DoubleArrayList probability, DoubleArrayList signal, short iteration) {
        short numEvents = (short) probability.size();
        double seprator = sum(signal)/2;
        DoubleArrayList positiveSignals = new DoubleArrayList();
        ShortArrayList positiveIndices = new ShortArrayList();
        DoubleArrayList negativeSignals = new DoubleArrayList();
        ShortArrayList negativeIndices = new ShortArrayList();
        
        
        short maxSignalIndex = maxIndex(signal);
        double w0 = 0.9;
        double w1 = 0.4;
        double w = (double) ((w0 - w1) * iteration * Math.sqrt(numEvents)) / 30;
        signal.set(maxSignalIndex, signal.get(maxSignalIndex) * w);
        
        
        for(short i = 0; i < numEvents; i++){
            if(signal.get(i) >= seprator) {
                positiveSignals.add(signal.get(i));
                positiveIndices.add(i);
            } else {
                negativeSignals.add(signal.get(i));
                negativeIndices.add(i);
            }
        }
        
        short positiveNum = (short) positiveSignals.size();
        double sumPositiveSignals = sum(positiveSignals);
        if(sumPositiveSignals > 0) {
            for(short i = 0; i < positiveNum; i++) {
                positiveSignals.set(i, (double) positiveSignals.get(i) / sumPositiveSignals);
            }
        }
        Collections.sort(positiveSignals);
        
        short negativeNum = (short) negativeSignals.size();
        double sumNegativeSignals = sum(negativeSignals);
        if(sumNegativeSignals > 0) {
            for(short i = 0; i < negativeNum; i++) {
                negativeSignals.set(i, (double) negativeSignals.get(i) / sumNegativeSignals);
            }
        } else {
            for(short i = 0; i < negativeNum; i++) {
                negativeSignals.set(i, (double) 1 / negativeNum);
            }
        }
        Collections.sort(negativeSignals);
        
       // System.out.println(negativeSignals);
        //System.out.println(negativeIndices);
        //System.out.println(positiveSignals);
        //System.out.println(positiveIndices);
        
        
        short index;
        for(short i = 0; i < negativeNum; i++) {
            index = negativeIndices.get(i);
            probability.set(index, (double) probability.get(index) * (1 - (negativeSignals.get(i) * beta)));
            
            for(short j = 0; j < numEvents; j++) {
                if(index != j) {
                    probability.set(j, (double) (((double) (negativeSignals.get(i) * beta)/(numEvents - 1)) + ((1 - (negativeSignals.get(i) * beta)) * probability.get(j))));
                }
            }                
//System.out.println(index + " " + probability + " " + sum(probability));            
        }
        
        
        
        for(short i = 0; i < positiveNum; i++) {
            index = positiveIndices.get(i);
            probability.set(index,(double) (probability.get(index) + ((positiveSignals.get(i) * alpha) * ( 1 - probability.get(index)))));
            for(short j = 0; j < numEvents; j++) {
                if(index != j) {
                    probability.set(j, probability.get(j) * (1 - (positiveSignals.get(i) * alpha)));
                }
            }         
        //System.out.println(index + " " + probability + " " + sum(probability));            
        }
        

        
        
        //System.out.println(positiveSignals);
        //System.out.println(sum(probability));
        //System.out.println(probability);
        
        //System.out.println(sum(probability) != 1);
        //double sumProbability = sum(probability);
        //if(sumProbability != 1) {
            //System.out.println("<---" + sum(probability) + " " +  (1 - sumProbability));
          //  probability.set(maxSignalIndex, (double) (probability.get(maxSignalIndex) + 1 - sumProbability));
            //System.out.println("--->" + sum(probability));
        //}
        //System.out.println(probability);
        
        /*
        for(short i = 0; i < negativeNum; i++) {
        
            beta_[kk] = sig_[kk] * algorithm->beta;
            nodes[j].probability[kk] = (1 - beta_[kk]) * nodes[j].probability[kk];
            for(k_ = 0; k_ < algorithm->num_partitions; k_++)
            {
                if(k_ != kk)
                    nodes[j].probability[k_] = (beta_[kk] / (algorithm->num_partitions - 1)) + (1 - beta_[kk]) * nodes[j].probability[k_];
            }
        }
        */
        
    }
    

    
    
    public static void main(String[] args)
    {
        int i = 0;
        int j = 0;
        
        int numberOfVertices = 5;
        int numberOfPartitions = 4;
        //System.out.println("numberOfPartitions " + numberOfPartitions);
        double alpha = 0.2;
        double beta = 0.02;
        short[] actions = new short[numberOfVertices];
        /*
        for(i = 0; i < numberOfVertices; i++)
        {
            actions[i] = (short) ThreadLocalRandom.current().nextInt(numberOfPartitions); //[0, max) inclusive 0 and exclusive max
            System.out.println("Action[" + i + "] = " + actions[i]);
        }
        
        
        double[][] probabilities = new double[numberOfVertices][numberOfPartitions]; 
        
        for(i = 0; i < numberOfVertices; i++)
        {
            Arrays.fill(probabilities[i], (double) 1/numberOfPartitions) ;
        }
        
        for(i = 0; i < numberOfVertices; i++)
        {
            System.out.print("Probabilities[" + i + "] = ");
            for(j = 0; j < numberOfPartitions; j++)
            {
                System.out.print(probabilities[i][j] + " ");
            }
            System.out.println();
        }
        
        
        short[] signals = new short[numberOfVertices];
        byte signalType = 2; // Binary \in [0, 1]
        for(i = 0; i < numberOfVertices; i++)
        {
            signals[i] = (short) ThreadLocalRandom.current().nextInt(signalType);
            System.out.println("Signals[" + i + "] = " + signals[i]);
        }
        
        for(i = 0; i < numberOfVertices; i++)
        {
            if(signals[i] == 0)
            {
                System.out.println("Reward");

                for(j = 0; j < numberOfPartitions; j++)
                {
                    if(j == actions[i])
                        probabilities[i][j] = probabilities[i][j] + (alpha * (1 - probabilities[i][j]));
                    else
                        probabilities[i][j] = (1 - alpha) * probabilities[i][j];
                }
            }
            else
            {
                System.out.println("Penalty");
                for(j = 0; j < numberOfPartitions; j++)
                {
                    if(j == actions[i])
                        probabilities[i][j] = (1 - beta) * probabilities[i][j];
                    else
                        probabilities[i][j] = (beta/(numberOfPartitions - 1)) + (1 - beta) * probabilities[i][j];
                }
            }
        }
        */
        /*
        for(i = 0; i < numberOfVertices; i++)
        {
            System.out.print("Probabilities[" + i + "] = ");
            for(j = 0; j < numberOfPartitions; j++)
            {
                System.out.print(probabilities[i][j] + " ");
            }
            System.out.println();
        }
        */
        
        List<Double> list = new ArrayList<Double>();
        short num = 8;
        for(i =0; i < num; i++) {
            list.add((double) (1 / num));
        }
        //System.out.println(list.get(0));
        
        DoubleArrayList LAProbability = new DoubleArrayList();
        for(i =0; i < num; i++) {
            LAProbability.add(0);
        }
        LAProbability.set(num - 3, .05);
        LAProbability.set(num - 2, .05);
       LAProbability.set(num - 1, .90);
            DoubleArrayList LASignal = new DoubleArrayList();
            for(i =0; i < num; i++) {
                LASignal.add(0);
                //LASignal1.add(10 * (double) Math.pow(i,2));
            }
            LASignal.set(num-2, 9000);
            LASignal.set(num-1, 3000);
       // LAProbability.set(0, 0);
        //LAProbability.set(1, 0);
       // LAProbability.set(2, 0);
       // LAProbability.set(3, 0);
       // LAProbability.set(4, 0);
        /*    
        LASignal1.add(100);
        System.out.println(LASignal1 + " " + maxIndex(LASignal1));    
        for(short k = 0; k < 30; k++) {             
            DoubleArrayList LASignal = new DoubleArrayList(LASignal1);
            System.out.println("+++++++++++++++++++++++++++++++++++++");
            System.out.println(LAProbability + " " + sum(LAProbability));
            System.out.println(actionSelection(LAProbability));
            weightedProbabilityUpdate(LAProbability, LASignal, k);
            System.out.println(LAProbability + " " + sum(LAProbability));
            System.out.println("---------------------------------------"); 
        }
        System.out.println(LASignal1 + " " + maxIndex(LASignal1));    
        */
        System.out.println(rouletteWheel(LAProbability));
        /*
        short k = 20;
        System.out.println(LASignal);
        System.out.println(LAProbability);
        weightedProbabilityUpdate(LAProbability, LASignal, k);
        System.out.println(LAProbability);
        
         ArrayList<Double> nfit = new ArrayList<Double>();
        nfit.add(1.0);
        nfit.add(3.0);
        nfit.add(5.0);
        nfit.add(21.0);
        nfit.add(5.0);
        nfit.add(45.0);
        System.out.println(nfit);

        ArrayList<Double> values = new ArrayList<Double>(nfit);
        System.out.println(values);
        values.set(1,1.1);
        System.out.println(nfit);
        System.out.println(values);
        int[] indices = new int[nfit.size()];
        indices[0] = 0;
        indices[1] = 1;
        indices[2] = 2;
        indices[3] = 3;
        indices[4] = 4;
        indices[4] = 5;

//int k = 0;
int l = 0;
double temp_value = 0.0;
int temp_index = 0;
int count = nfit.size();
for(k = 0; k < count; k++) {
    for(l = k + 1; l < count; l++) {
        if(values.get(k) > values.get(l)) {
                temp_value = values.get(k);
                values.set(k,values.get(l));
                values.set(l,temp_value);
                temp_index = indices[k];
                indices[k] = indices[l];
                indices[l] = temp_index;
}
}
}

System.out.println(values);
for(i = 0; i < count; i++)
   System.out.print(indices[i]);
   System.out.println();     
        */
        
        
        
        
        System.out.println("++++++++++++++++++++++++++++++++++");
        
        
        
    }

}
