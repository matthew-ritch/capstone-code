//much of source is adapted from https://github.com/ejoebstl/Android-Sensor-Log/blob/master/app/src/main/java/io/iam360/sensorlog/MainActivity.java

package capstone.wheelchairmonitor;


import android.annotation.SuppressLint;
import android.content.Context;
import android.content.Intent;
import android.content.pm.PackageInfo;
import android.content.pm.PackageManager;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener2;
import android.hardware.SensorManager;
import android.os.Environment;
import android.os.Bundle;
import android.os.SystemClock;
import android.util.Log;
import android.view.KeyEvent;
import android.view.MotionEvent;
import android.view.View;
import android.view.inputmethod.EditorInfo;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

import androidx.appcompat.app.AppCompatActivity;

import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.matrix.dense.Basic2DMatrix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.text.Format;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;


public class MainActivity extends AppCompatActivity implements SensorEventListener2 {

    SensorManager manager;
    Button buttonCalib;
    Button buttonRecord;
    Button buttonProcess;
    Button buttonSelect;
    Button buttonStop;
    EditText editText;
    TextView resultsTextView;
    static File toRead = null;
    static Boolean proceed = false;
    String toReadFold;
    String toReadFile;
    int Fs = 100; // 1/s
    int t_calib = 3; // s
    boolean isRunning;
    boolean magFirst;
    final String TAG = "SensorLog";
    FileWriter writer;
    double[] calib_acc = {0, 0, 0};
    double[] calib_mag = {0, 0, 0};
    static long stoptime = System.currentTimeMillis();
    int PICKFILE_RESULT_CODE = 100;
    static String runNam = "";


    @SuppressLint("ClickableViewAccessibility")
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        isRunning = false;
        magFirst = false;
        manager = (SensorManager) getSystemService(Context.SENSOR_SERVICE);
        buttonCalib = (Button)findViewById(R.id.button1);
        buttonRecord = (Button)findViewById(R.id.button2);
        buttonSelect = (Button)findViewById(R.id.button3);
        buttonStop = (Button)findViewById(R.id.buttonStop);
        buttonProcess = (Button)findViewById(R.id.button_process);
        editText = (EditText) findViewById(R.id.editText);
        resultsTextView = (TextView) findViewById(R.id.textView);
        //initial state
        buttonProcess.setEnabled(false);

        buttonCalib.setOnTouchListener(new View.OnTouchListener() {
            @Override
            public boolean onTouch(View view, MotionEvent motionEvent) {
                if (motionEvent.getAction() == MotionEvent.ACTION_DOWN) {
                    buttonCalib.setEnabled(false);
                    buttonRecord.setEnabled(false);
                    buttonProcess.setEnabled(false);

                    Log.d(TAG, "Loading the writer...\n");
                    try {

                        writer = new FileWriter(new File(getStorageDir(), "sensors_" + "calibration" + ".csv"));
                    } catch (IOException e) {
                        e.printStackTrace();
                        return false;
                    }
                    try {
                        manager.registerListener(MainActivity.this, manager.getDefaultSensor(Sensor.TYPE_ACCELEROMETER), (10 ^ 6) / Fs);
                        manager.registerListener(MainActivity.this, manager.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD), (10 ^ 6) / Fs);
                        Log.d(TAG, "manager registered...\n");
                    } catch (Exception e) {
                        return false;
                    }
                    Log.d(TAG, "Writing calibration data to " + getStorageDir());
                    // Store the start time
                    stoptime=System.currentTimeMillis()+t_calib*1000;
                    isRunning = true;
                    return true;
                }
                else {
                    return true;
                }
            }



        });

        buttonRecord.setOnTouchListener(new View.OnTouchListener() {
            @Override
            public boolean onTouch(View view, MotionEvent motionEvent) {
                if (motionEvent.getAction() == MotionEvent.ACTION_DOWN) {
                    buttonCalib.setEnabled(false);
                    buttonRecord.setEnabled(false);
                    buttonProcess.setEnabled(false);

                    Log.d(TAG, "Loading the writer...\n");
                    try {
                        runNam=editText.getText().toString();
                        writer = new FileWriter(new File(getStorageDir(), "sensors_recording_" + runNam + ".csv"));
                    } catch (IOException e) {
                        e.printStackTrace();
                        return false;
                    }
                    try {
                        manager.registerListener(MainActivity.this, manager.getDefaultSensor(Sensor.TYPE_ACCELEROMETER), (10 ^ 6) / Fs);
                        manager.registerListener(MainActivity.this, manager.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD), (10 ^ 6) / Fs);
                        Log.d(TAG, "manager registered...\n");
                    } catch (Exception e) {
                        return false;
                    }
                    Log.d(TAG, "Writing calibration data to " + getStorageDir());
                    // Store the start time
                    isRunning = true;
                    return true;
                }
                else {
                    return true;
                }
            }



        });

        buttonStop.setOnTouchListener(new View.OnTouchListener() {
            @Override
            public boolean onTouch(View view, MotionEvent motionEvent) {
                if (motionEvent.getAction() != MotionEvent.ACTION_DOWN) {
                    return false;
                }
                isRunning = false;
                try {
                    writer.close();//close writer
                } catch (IOException e) {
                    e.printStackTrace();
                }
                buttonCalib.setEnabled(true);
                buttonRecord.setEnabled(true);
                buttonProcess.setEnabled(true);
                manager.flush(MainActivity.this);
                manager.unregisterListener(MainActivity.this);
                // check calibration
                calib_acc = getCalib_acc();
                Log.d("check calibration acc", String.valueOf(calib_acc[0]) + " " + String.valueOf(calib_acc[1]) + " " + String.valueOf(calib_acc[2]));
                calib_mag = getCalib_mag();
                Log.d("check calibration mag", String.valueOf(calib_mag[0]) + " " + String.valueOf(calib_mag[1]) + " " + String.valueOf(calib_mag[2]));
                Log.d(TAG, "previous task stopped\n");
                return true;
            }
        });

        buttonSelect.setOnTouchListener(new View.OnTouchListener() {
            @Override
            public boolean onTouch(View view, MotionEvent motionEvent) {
                if (motionEvent.getAction() != MotionEvent.ACTION_DOWN) {
                    return false;
                }
                buttonProcess.setEnabled(true);
                /*//load in calibration file
                calib_acc = getCalib_acc();
                Log.d("check calibration", String.valueOf(calib_acc[0]) + " " + String.valueOf(calib_acc[1]) + " " + String.valueOf(calib_acc[2]));
                calib_mag = getCalib_mag();
                Log.d("check calibration", String.valueOf(calib_mag[0]) + " " + String.valueOf(calib_mag[1]) + " " + String.valueOf(calib_mag[2]));*/

                //select file
                Intent intent = new Intent(Intent.ACTION_GET_CONTENT);
                intent.setType("*/*");
                startActivityForResult(intent, PICKFILE_RESULT_CODE); //this changes toRead
                Log.d("check","made it past safr");
                proceed = false;//reset proceed
                return true;
            }
        });

        buttonProcess.setOnTouchListener(new View.OnTouchListener() {
            @Override
            public boolean onTouch(View view, MotionEvent motionEvent) {
                if (motionEvent.getAction() != MotionEvent.ACTION_DOWN) {
                    return false;
                }
                //load in calibration file
                calib_acc = getCalib_acc();
                Log.d("check calibration", String.valueOf(calib_acc[0]) + " " + String.valueOf(calib_acc[1]) + " " + String.valueOf(calib_acc[2]));
                calib_mag = getCalib_mag();
                Log.d("check calibration", String.valueOf(calib_mag[0]) + " " + String.valueOf(calib_mag[1]) + " " + String.valueOf(calib_mag[2]));
                // initialize matrices
                Basic2DMatrix accMat = Basic2DMatrix.zero(1,3);//0,3
                Basic2DMatrix magMat = Basic2DMatrix.zero(1,3);//0,3
                //read in file
                //http://la4j.org/apidocs/org/la4j/Matrix.html#insertRow(int,%20org.la4j.Vector)
                List<Double> T = new ArrayList<>();
                if (toRead!=null) {
                    try {
                        accMat = Basic2DMatrix.zero(1,3);//0,3
                        magMat = Basic2DMatrix.zero(1,3);//0,3
                        //Log.d("processing", "to read in not null");
                        BufferedReader br = new BufferedReader(new FileReader(toRead));
                        String line;
                        int i_acc=0; int i_mag=0;
                        while ((line = br.readLine()) != null) {
                            //Log.d("processing line", line);
                            String[] values = line.split(",");
                            //extract vectors of measurements


                            if (values[1].equals("ACC")) {
                                accMat= (Basic2DMatrix) accMat.insertRow(i_acc, Vector.fromArray(new double[]{Double.parseDouble(values[2]), Double.parseDouble(values[3]), Double.parseDouble(values[4])}));
                                //Log.d("nAcc prog", String.valueOf(accMat.rows()));
                                T.add(Double.parseDouble(values[0]));
                                i_acc+=1;
                                continue;
                            }
                            if (values[1].equals("MAG")) {
                                magMat= (Basic2DMatrix) magMat.insertRow(i_mag, Vector.fromArray(new double[]{Double.parseDouble(values[2]), Double.parseDouble(values[3]), Double.parseDouble(values[4])}));
                                //Log.d("nMag prog", String.valueOf(magMat.rows()));
                                i_mag+=1;
                                continue;
                            }

                        }
                        br.close();
                    } catch (IOException e) {
                        //could not read file
                    }
                    //do processing.

                    //find effective sampling rates - acc is 500 hz and mag is 100 hz..
                    int nAcc=accMat.rows();
                    Log.d("nAcc post", String.valueOf(nAcc));
                    int nMag=magMat.rows();
                    Log.d("nMag post", String.valueOf(nMag));
                    double len = (T.get(nAcc - 2) - T.get(0))/Math.pow(10,9);//minus two here because of first row of zeros on acc and mag
                    Log.d("len", String.valueOf(len));
                    double Fs_acc = nAcc / len;
                    Log.d("Fs_acc", String.valueOf(Fs_acc));
                    double Fs_mag = nMag / len;
                    Log.d("Fs_mag", String.valueOf(Fs_mag));

                    int sampratio = Math.round(nAcc/nMag);
                    //downsample - done in naive way. very consistent 5 mag to 1 acc. la4j removeRow
                    int[] slicing = new int[(int) Math.floor(nAcc/sampratio)];
                    //List<Integer> slicing = new ArrayList<>();
                    for (int k = 0; k++ < Math.floor(nAcc/sampratio)-1;) {slicing[k]=(k*sampratio);}
                    //accMat = accMat;
                    accMat = (Basic2DMatrix) accMat.select(slicing, new int[]{0, 1, 2});


                    //make lengths work together - no evidence for failure here yet!
                    while (accMat.rows()<magMat.rows()) {magMat.removeLastRow();}
                    while (accMat.rows()>magMat.rows()) {accMat.removeLastRow();}
                    nMag=magMat.rows();
                    Log.d("nMag down", String.valueOf(nMag));
                    nAcc=accMat.rows();
                    Log.d("nAcc down", String.valueOf(nAcc));

                    //decide which time vector to use, mag or acc. mag is less noisy (?). will make own idealized version
                    double Fs_effective = nAcc/len;
                    //tvec= //for loop

                    //remove gravity from acceleration
                    Basic2DMatrix rotG = new Basic2DMatrix(accMat.rows(),3);
                    //use calib_acc, calib_mag
                    Vector magorient = Vector.fromArray(calib_mag);
                    magorient.divide(magorient.norm());
                    Vector g = Vector.fromArray(calib_acc);
                    //step through and subtract rotated gravity from each row
                    for (int i=0; i++<accMat.rows();){
                        Vector a = magMat.getRow(i).divide(magMat.getRow(i).norm());
                        Matrix rotator = RU(a, magorient);
                        accMat.setRow(i, accMat.getRow(i).subtract(rotator.multiply(g)) );
                    }

                    Vector sum_params = magMat.getColumn(0).add(magMat.getColumn(1).add(magMat.getColumn(2))).add(accMat.getColumn(0).add(accMat.getColumn(1).add(accMat.getColumn(2))));
                    Matrix all_params = accMat.insertColumn(3, magMat.getColumn(0)).insertColumn(4, magMat.getColumn(1)).insertColumn(5, magMat.getColumn(2));

                    FFT fft = new FFT(256);
                    double [] x = new double[256];
                    double [] y = new double[256];
                    double [] p = new double[128];
                    double f_dom = 0;
                    double p_dom = 0;

                    int n_sampleSets = (int) Math.floor(accMat.rows() / 256);
                    //Vector weights = Vector.fromArray(new double[] {});
                    Vector simpweights = Vector.fromArray(new double[]{-7.927104e-03, -6.632750e-02, -1.597056e-02, -1.376609e-01, -1.548602e-01, -2.369039e-01, 3.380047e-03, -7.083327e-05, 4.825574e-03, 1.869104e-01, 1.792355e-01, -1.785958e-01, 1.239226e-02, -1.799500e-02, -7.220091e-02, -1.246121e-01, -4.310982e-01, -3.226717e-01, 3.999075e-02, 1.961273e-02, -5.612607e-02, -2.697329e-01, 2.231319e+00, 2.838654e-01, -7.530079e-04, 2.698552e-04, 8.214414e-04, -7.517987e-02, -7.391309e-03, 8.287098e-02, -8.545501e-02, 1.433272e-01, 4.602780e-02, 1.733912e+00, -6.382943e+00, -2.634976e+00, 7.760145e-01,});
                    Vector actweights = Vector.fromArray(new double[]{5.908857e-02,1.249869e-01,-5.220351e-02,9.694006e-02,-6.417038e-01,1.635168e+00,3.712133e-02,3.700235e-04,4.273163e-02,1.129824e+00,1.463459e-01,7.889935e-01,-6.956490e-02,-8.952888e-02,7.942882e-02,5.986386e-01,-3.906612e-01,-4.999971e-01,4.694541e-01,1.650772e-01,-3.449564e-01,-1.777032e+00,-1.017185e+01,1.119613e+00,-8.652106e-03,2.023478e-03,2.612582e-03,2.647245e-01,-3.258171e-01,3.745326e-01,-6.471418e-01,-5.648614e-01,3.513456e-01,6.502583e+00,2.027482e+01,-5.942074e+00,4.365293e+00});
                    double active_n = 0 ;
                    double[] types_active = new double[3];
                    for (int i = 0; i++<n_sampleSets;){
                        //slice(int fromRow, int fromColumn, int untilRow, int untilColumn)
                        Matrix samples = all_params.slice(i*256,0,(i+1)*256,6);
                        Vector sum_samples = sum_params.slice(i*256, (i+1)*256);
                        for (int j=0;j++<256;){x[j] = sum_samples.get(j);}
                        //do fft
                        fft.fft(x, y);
                        double f_dom_this = 0;
                        for (int j=0;j++<128;){
                            double p_this = Math.pow(x[j],2) + Math.pow(y[j],2);
                            if (p_this > p_dom){
                                p_dom = p_this;
                                f_dom_this = j*(Fs_effective/256);
                            }
                        }


                        //do regression stuff
                        Vector features = getFeatures(samples, Fs_effective);
                        //calculate active vs inactive classification
                        double active_stat = features.hadamardProduct(simpweights).sum();
                        if (active_stat>1.5){
                            active_n+=1;
                            f_dom+=f_dom_this;
                            //calculate activity classification
                            double type_active = features.hadamardProduct(actweights).sum();
                            if (type_active<1.5){
                                types_active[0]+=1;
                            }
                            else if (type_active>2.5){
                                types_active[2]+=1;
                            }
                            else {
                                types_active[1]+=1;
                            }

                        }

                    }
                    //display results
                    int pct_active = (int) Math.round(100*active_n/n_sampleSets);
                    int [] pcts_active_types = new int[] {(int) Math.round(100 * types_active[0] / active_n),(int) Math.round(100 * types_active[1] / active_n),(int) Math.round(100 * types_active[2] / active_n)};
                    f_dom=f_dom/active_n;
                    String stringResults = String.format("File processed.\nThe user was %i%% active and %i%% inactive during this recorded session.\n" +
                            "While actively pushing, the user used proper form %i%% of the time, was lifting %i%% of the time, was choppy %i%% of the time, and averaged a frequency of %d Hz.",
                            pct_active, 100-pct_active, pcts_active_types[1],pcts_active_types[0],pcts_active_types[2],f_dom);
                    resultsTextView.setText(stringResults);


                }



                return true;
            }
        });






    }

    Vector getFeatures(Matrix samples, double fs){
        //mean, stdev, rms,
        int n = samples.rows();
        //Vector features = Vector.fromArray(new double[] {samples.getColumn(0).sum()/n, samples.getColumn(1).sum()/n, samples.getColumn(2).sum()/n,samples.getColumn(3).sum()/n, samples.getColumn(4).sum()/n,samples.getColumn(5).sum()/n,samples.getColumn(0).subtract(means.get(0)).sum()/n, samples.getColumn(1).subtract(means.get(1)).sum()/n, samples.getColumn(2).subtract(means.get(2)).sum()/n,samples.getColumn(3).subtract(means.get(3)).sum()/n, samples.getColumn(4).subtract(means.get(4)).sum()/n,samples.getColumn(5).subtract(means.get(5)).sum()/n});
        //Vector stds = Vector.fromArray(new double[] {samples.getColumn(0).subtract(means.get(0)).sum()/n, samples.getColumn(1).subtract(means.get(1)).sum()/n, samples.getColumn(2).subtract(means.get(2)).sum()/n,samples.getColumn(3).subtract(means.get(3)).sum()/n, samples.getColumn(4).subtract(means.get(4)).sum()/n,samples.getColumn(5).subtract(means.get(5)).sum()/n});
        Matrix ints = integrateCols(samples).multiply(1/fs);
        Matrix derivs = integrateCols(samples).multiply(fs);
        samples = samples.insert(derivs, 0, 6, 256, 6).insert(ints, 0, 6, 256, 6);

        Vector means = sumCols(samples).divide(n);
        Vector stds = sumCols(subtractCols(samples, means)).divide(n);
        double[] features = new double [37];
        features[36]=1;
        for (int j=0;j++<18;){features[j] = means.get(j);}
        for (int j=0;j++<18;){features[j+18] = stds.get(j);}
        //TODO add cross correlations

        return Vector.fromArray(features);
    }

    Vector sumCols(Matrix mat){
        double outs [] = new double [mat.columns()];
        for (int i=0;i++<mat.columns();){
            outs[i] = mat.getColumn(i).sum();
        }
        return Vector.fromArray(outs);
    }

    Matrix subtractCols(Matrix mat, Vector v){
        for (int i=0;i++<mat.rows();){
            mat.setRow(i,mat.getRow(i).subtract(v));
        }
        return mat;
    }



    Matrix integrateCols(Matrix mat){
        Matrix out = mat.copy();
        Vector runSum = Vector.zero(mat.columns());
        for (int i = 0; i++< mat.rows();){
            runSum = runSum.add(mat.getRow(i));
            out.setRow(i, runSum);
        }
        return out;
    }

    Matrix differentiateCols(Matrix mat){
        Matrix out = Matrix.zero(mat.rows(),mat.columns());
        for (int i = 1; i++< mat.rows();){
            out.setRow(i, mat.getRow(i).subtract(mat.getRow(i-1)));
        }
        return out;
    }

    Matrix ssc(Vector v){
        // ssc=@(x)[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
        double[][] arr = {{0, -v.get(2), v.get(1)},{v.get(2), 0, -v.get(0)},{-v.get(1), v.get(0), 0}};
        return Matrix.from2DArray(arr);
    }

    Matrix RU(Vector A, Vector B){
        //RU = @(A,B) eye(3) + ssc(cross(A,B)) + ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2);
        Matrix added = Basic2DMatrix.identity(3).add(  ssc(cross(A, B))  .  add( ssc(cross(A, B)).power(2).multiply((1 - A.hadamardProduct(B).sum()) / (cross(A, B).norm())))  );
        return added;
    }

    Vector cross(Vector A, Vector B){
    //⟨a1b2−a2b1,a2b0−a0b2,a0b1−a1b0⟩
        return Vector.fromArray(new double[]{A.get(1) * B.get(2) - A.get(2) * B.get(1), A.get(2) * B.get(0) - A.get(0) * B.get(2), A.get(0) * B.get(1) - A.get(1) * B.get(0)});
    }


    protected void onActivityResult(int requestCode, int resultCode, Intent data) {
        if (requestCode == PICKFILE_RESULT_CODE) {

            String FilePath = data.getData().getPath();
            Log.d("FilePath", FilePath);
            String toReadFile = data.getData().getLastPathSegment();
            Log.d("toReadFile", toReadFile);
            int lastPos = FilePath.length() - toReadFile.length();
            String toReadFold = FilePath.substring(0, lastPos);
            Log.d("toReadFold", toReadFold);

            int ind = FilePath.lastIndexOf("/");
            String fileName = FilePath.substring(ind+1);//start one past for file name
            Log.d("fileName", fileName);


            
            toRead = new File(getStorageDir(), fileName);
            proceed = true;
            return;
        }
        super.onActivityResult(requestCode, resultCode, data);
    }

    private String getStorageDir() {
        return this.getExternalFilesDir(null).getAbsolutePath();
        //  return "/storage/emulated/0/Android/data/com.iam360.sensorlog/";
    }

    public double [] getCalib_acc() {
        double[] acc = {0, 0, 0};
        //Get the text file
        File readfile = new File(getStorageDir(), "sensors_" + "calibration" + ".csv");
        //Read text from file
        int i = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(readfile));
            String line;

            while ((line = br.readLine()) != null) {
                //Log.d("getCalib line", line);
                String[] values = line.split(",");
                //make sure line is an acceleration line
                //Log.d("getCalib label", values[1]);
                if  (!values[1].equals("ACC")) {
                    continue;
                }
                i+=1; //record how many values have been added
                //Log.d("getCalib", String.valueOf(i));
                acc[0]+=Float.parseFloat(values[2]);
                acc[1]+=Float.parseFloat(values[3]);
                acc[2]+=Float.parseFloat(values[4]);

            }
            br.close();
        }
        catch (IOException e) {
            //could not read file
        }
        double[] calib_acc = {acc[0]/i, acc[1]/i, acc[2]/i};
        return calib_acc;
    }

    public double [] getCalib_mag() {
        double[] mag = {0, 0, 0};
        //Get the text file
        File readfile = new File(getStorageDir(), "sensors_" + "calibration" + ".csv");
        //Read text from file
        int i = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(readfile));
            String line;

            while ((line = br.readLine()) != null) {

                String[] values = line.split(",");
                //make sure line is an acceleration line
                if  (!values[1].equals("MAG")) {
                    continue;
                }
                i+=1; //record how many values have been added
                mag[0]+=Float.parseFloat(values[2]);
                mag[1]+=Float.parseFloat(values[3]);
                mag[2]+=Float.parseFloat(values[4]);

            }
            br.close();
        }
        catch (IOException e) {
            //could not read file
        }
        double[] calib_mag = {mag[0]/i, mag[1]/i, mag[2]/i};
        return calib_mag;

    }

    @Override
    public void onFlushCompleted(Sensor sensor) {

    }

    //TODO make so runs for a specified amount of time
    @Override
    public void onSensorChanged(SensorEvent evt) {
        if(isRunning ) {
            //Log.d("sensorwriting", "passed isRunning\n");
            try {
                switch(evt.sensor.getType()) {
                    case Sensor.TYPE_ACCELEROMETER:
                        if (!magFirst){break;} //mag takes longer to spin up
                        writer.write(String.format("%d,ACC,%f,%f,%f\n", evt.timestamp, evt.values[0], evt.values[1], evt.values[2]));
                        //Log.d("sensorwriting", "wrote acc\n");
                        break;
                    case Sensor.TYPE_MAGNETIC_FIELD:
                        if (!magFirst){magFirst = true;} //mag takes longer to spin up
                        //timestamp in nanoseconds
                        writer.write(String.format("%d,MAG,%f,%f,%f\n", evt.timestamp, evt.values[0], evt.values[1], evt.values[2]));
                        //Log.d("sensorwriting", "wrote MAG\n");
                        break;
                }
            } catch (IOException e) {
                e.printStackTrace();
                Log.d(TAG, "writer failed\n");
            }
        }
        else {
            buttonStop.performClick();
        }
    }

    @Override
    public void onAccuracyChanged(Sensor sensor, int i) {

    }
}