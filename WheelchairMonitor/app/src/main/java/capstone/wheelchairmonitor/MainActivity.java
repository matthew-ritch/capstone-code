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
                    buttonSelect.setEnabled(false);

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
                    buttonSelect.setEnabled(false);

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
                buttonSelect.setEnabled(true);
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

                buttonCalib.setEnabled(false);
                buttonRecord.setEnabled(false);
                buttonProcess.setEnabled(false);
                buttonSelect.setEnabled(false);
                resultsTextView.setText("Processing, please wait...");

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
                        int i_acc=0; int i_mag=0; int ii=0;
                        while ((line = br.readLine()) != null) {
                            //Log.d("processing line", line);
                            String[] values = line.split(",");
                            //extract vectors of measurements
                            //drop first 100 samples
                            if (ii<100){
                                ii++;
                                continue;
                            }



                            if (values[1].equals("ACC")) {
                                //if (i_acc%5 == 0){
                                    accMat= (Basic2DMatrix) accMat.insertRow(i_acc, Vector.fromArray(new double[]{Double.parseDouble(values[2]), Double.parseDouble(values[3]), Double.parseDouble(values[4])}));
                                    //Log.d("nAcc prog", String.valueOf(accMat.rows()));
                                    T.add(Double.parseDouble(values[0]));
                                //}
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
                    //Log.d("nAcc post", String.valueOf(nAcc));
                    int nMag=magMat.rows();
                    //Log.d("nMag post", String.valueOf(nMag));
                    double len = (T.get(nAcc - 2) - T.get(0))/Math.pow(10,9);//minus two here because of first row of zeros on acc and mag
                    //Log.d("len", String.valueOf(len));
                    double Fs_acc = nAcc / len;
                    //Log.d("Fs_acc", String.valueOf(Fs_acc));
                    double Fs_mag = nMag / len;
                    //Log.d("Fs_mag", String.valueOf(Fs_mag));

                    int sampratio = 5;


                    //downsample - done in naive way. very consistent 5 mag to 1 acc. la4j removeRow
                    int[] slicing = new int[(int) Math.floor(nAcc/sampratio)];
                    //Log.d("slicing", "make slicing array done.");
                    //List<Integer> slicing = new ArrayList<>();
                    for (int k = 0; k < Math.floor(nAcc/sampratio); k++) {slicing[k]=(k*sampratio);}
                    //accMat = accMat;
                    accMat = (Basic2DMatrix) accMat.select(slicing, new int[]{0, 1, 2});
                    //
                    Log.d("slicing", "slicing done.");

                    //Mag=magMat.rows();
                    //Log.d("nMag down", String.valueOf(nMag));
                    //nAcc=accMat.rows();
                    //Log.d("nAcc down", String.valueOf(nAcc));


                    //make lengths work together - no evidence for failure here yet!
                    while (accMat.rows()<magMat.rows()) {magMat = (Basic2DMatrix) magMat.removeLastRow();Log.d("removing", "removed one from magMat");}
                    while (accMat.rows()>magMat.rows()) {accMat = (Basic2DMatrix) accMat.removeLastRow();Log.d("removing", "removed one from accMat");}
                    nMag=magMat.rows();
                    Log.d("nMag downsampled", String.valueOf(nMag));
                    nAcc=accMat.rows();
                    Log.d("nAcc downsampled", String.valueOf(nAcc));

                    //decide which time vector to use, mag or acc. mag is less noisy (?). will make own idealized version
                    double Fs_effective = nAcc/len;
                    //tvec= //for loop

                    //remove gravity from acceleration
                    //Basic2DMatrix rotG = new Basic2DMatrix(accMat.rows(),3);
                    //use calib_acc, calib_mag
                    Vector magorient = Vector.fromArray(calib_mag);
                    magorient=magorient.divide(magorient.norm());
                    Vector g = Vector.fromArray(calib_acc);
                    //step through and subtract rotated gravity from each row
                    //fixed: why does this result in diverging acceleration values past g? not just random rotation error (see accmatx    346 going past 20 in some test cases)
                    if (true) {
                        for (int i = 0; i < accMat.rows(); i++) {
                            Vector a = magMat.getRow(i).divide(magMat.getRow(i).norm());
                            magMat.setRow(i, a);
                            if (false) {
                                Matrix rotator = RU(a, magorient);
                                //Log.d("rotator", rotator.toString());
                                accMat.setRow(i, accMat.getRow(i).subtract(rotator.multiply(g)));
                            }
                        }
                    }
                    //Log.d("magMat 344", magMat.toString());
                    //Log.d("accMat 346", accMat.toString());

                    Vector sum_params = magMat.getColumn(0).add(magMat.getColumn(1).add(magMat.getColumn(2))).add(accMat.getColumn(0).add(accMat.getColumn(1).add(accMat.getColumn(2))));
                    Basic2DMatrix all_params = new Basic2DMatrix(accMat.rows(), 6);
                    //combine matrices
                    for (int i = 0; i<3; i++){all_params.setColumn(i,accMat.getColumn(i));}
                    for (int i = 0; i<3; i++){all_params.setColumn(i+3,magMat.getColumn(i));}

                    FFT fft = new FFT(256);
                    double [] x = new double[256];
                    double [] y = new double[256];
                    //double [] p = new double[128];
                    double f_dom = 0;
                    double p_dom = 0;

                    int n_sampleSets = (int) Math.floor(accMat.rows() / 256);
                    //Vector weights = Vector.fromArray(new double[] {});
                    Vector simpweights = Vector.fromArray(new double[]{5.755775E-03, -6.699398E-02, -1.445151E-02, -3.343928E-02, 2.393271E-01, 1.270501E-02, 2.398644E-03, -4.778000E-03, 2.975728E-03, 1.288486E-01, -2.109809E-01, 2.745229E-01, 7.096885E-03, -5.360486E-02, 5.685067E-03, 1.661659E-01, -8.727464E-01, 3.387613E-01, -7.854814E-03, 4.691318E-02, -7.993375E-02, 2.528162E-01, 2.816685E+00, -6.074269E-01, 3.497633E-04, -7.123845E-04, 7.866524E-04, -1.269227E-01, -3.467171E-04, 1.098068E-01, -1.683645E-02, 2.059982E-01, 1.041541E-01, 3.357326E-01, -8.450184E+00, 2.392375E-01, 1.098618E+00});
                    Vector actweights = Vector.fromArray(new double[]{2.047621E-01, 2.007387E-01, -9.797865E-02, 5.776776E-01, -1.580921E+00, 5.020793E-01, 2.603543E-03, -2.583976E-05, 1.564483E-03, 1.829132E-01, -1.886063E-02, 3.414378E-01, -6.893037E-02, 2.165704E-02, -6.320788E-02, -1.652131E-02, -1.048148E+00, 7.957762E-02, 1.541733E-01, 2.827244E-01, -3.819768E-01, -8.834009E-01, -6.253969E+00, -7.779593E-01, -2.824386E-03, -1.561134E-03, 2.726407E-03, 3.770149E-01, -6.818273E-01, 2.303266E-01, -1.017210E-02, -7.628410E-01, 1.733548E-01, 4.624819E+00, 5.340712E+00, 8.938623E-01, 5.877527E+00});
                    Log.d("simpweights",simpweights.toString());
                    double active_n = 0 ;
                    double[] types_active = new double[3];
                    for (int i = 0; i<n_sampleSets;i++){
                        Log.d("samplesets", "completed sample " + String.valueOf(i+1)+" of total "+ String.valueOf(n_sampleSets));
                        //slice(int fromRow, int fromColumn, int untilRow, int untilColumn)
                        Matrix samples = all_params.slice(i*256,0,(i+1)*256,6);
                        Vector sum_samples = sum_params.slice(i*256, (i+1)*256);
                        for (int j=0;j<256;j++){x[j] = sum_samples.get(j);}
                        //do fft
                        fft.fft(x, y);
                        double f_dom_this = 0;
                        for (int j=0;j<26;j++){
                            double p_this = Math.pow(x[j],2) + Math.pow(y[j],2);
                            if (p_this > p_dom){
                                p_dom = p_this;
                                f_dom_this = j*(Fs_effective/256);
                            }
                        }


                        //do regression stuff
                        Vector features = getFeatures(samples, Fs_effective);
                        //Log.d("features", features.toString());
                        //calculate active vs inactive classification
                        double active_stat = features.copy().hadamardProduct(simpweights).sum();
                        //Log.d("features", features.toString());
                        Log.d("active_stat", String.valueOf((active_stat)));
                        //calculate activity classification
                        double type_active = features.copy().hadamardProduct(actweights).sum();
                        //Log.d("features", features.toString());
                        Log.d("type_active", String.valueOf(type_active));

                        if (active_stat>1.5){
                            active_n+=1;
                            f_dom+=f_dom_this;

                            //TODO move fft to this part
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
                    int pct_active = (int) Math.round(100.00*active_n/n_sampleSets);
                    int [] pcts_active_types = new int[] {(int) Math.round(100 * types_active[0] / active_n),(int) Math.round(100.00 * types_active[1] / active_n),(int) Math.round(100.00 * types_active[2] / active_n)};
                    if (active_n<1){active_n=1;}
                    f_dom= f_dom/active_n;
                    //TODO add all rest case to this
                    String stringResults;
                    if (pct_active>0) {
                        stringResults = String.format("File processed.\nThe user was %d%% active and %d%% inactive during this recorded session.\n" +
                                        "While actively pushing, the user used proper form %d%% of the time, was lifting %d%% of the time, was choppy %d%% of the time, and averaged a frequency of %f Hz.",
                                pct_active, 100 - pct_active, pcts_active_types[1], pcts_active_types[0], pcts_active_types[2], f_dom);
                    }
                    else {
                        stringResults = String.format("File processed.\nThe user was not actively pushing during this recorded session.\n");
                    }
                    resultsTextView.setText(stringResults);
                    buttonCalib.setEnabled(true);
                    buttonRecord.setEnabled(true);
                    buttonProcess.setEnabled(true);
                    buttonSelect.setEnabled(true);
                }
                return true;
            }
        });


    }
    Vector getFeatures(Matrix samples, double fs){
        fs = 100.000;
        //mean, stdev, rms,
        Matrix matIn = samples.copy();
        double n = samples.rows();
        //Log.d("samples 440", samples.toString());
        Vector means = sumCols(samples).divide(n);
        //Log.d("means 441", means.toString());
        Basic2DMatrix mat_centered = (Basic2DMatrix) subtractCols(samples, means);
        //Vector features = Vector.fromArray(new double[] {samples.getColumn(0).sum()/n, samples.getColumn(1).sum()/n, samples.getColumn(2).sum()/n,samples.getColumn(3).sum()/n, samples.getColumn(4).sum()/n,samples.getColumn(5).sum()/n,samples.getColumn(0).subtract(means.get(0)).sum()/n, samples.getColumn(1).subtract(means.get(1)).sum()/n, samples.getColumn(2).subtract(means.get(2)).sum()/n,samples.getColumn(3).subtract(means.get(3)).sum()/n, samples.getColumn(4).subtract(means.get(4)).sum()/n,samples.getColumn(5).subtract(means.get(5)).sum()/n});
        //Vector stds = Vector.fromArray(new double[] {samples.getColumn(0).subtract(means.get(0)).sum()/n, samples.getColumn(1).subtract(means.get(1)).sum()/n, samples.getColumn(2).subtract(means.get(2)).sum()/n,samples.getColumn(3).subtract(means.get(3)).sum()/n, samples.getColumn(4).subtract(means.get(4)).sum()/n,samples.getColumn(5).subtract(means.get(5)).sum()/n});

        Basic2DMatrix ints = (Basic2DMatrix) integrateCols(mat_centered).divide(fs);
        //Log.d("ints",ints.toString());
        Matrix derivs = differentiateCols(mat_centered).multiply(fs);
        Matrix trip = new Basic2DMatrix(256,18);
        trip = trip.insert(matIn,0,0,256,6).insert(derivs, 0, 6, 256, 6).insert(ints, 0, 12, 256, 6);
        Log.d("trip 451", trip.toString());
        //Log.d("samples 451", samples.toString());
        means = sumCols(trip).divide(n);
        //Log.d("means 450", means.toString());
        Basic2DMatrix diffs = (Basic2DMatrix) subtractCols(trip, means);



        Vector stds = sumCols(  diffs.hadamardProduct(diffs)  ).divide(n-1);
        Log.d("stds 456", stds.toString());
        double[] features = new double [37];
        features[36]=1;
        for (int j=0;j<18;j++){features[j] = means.get(j);}
        for (int j=0;j<18;j++){features[j+18] = Math.sqrt(stds.get(j));}
        //TODO add cross correlations

        return Vector.fromArray(features);
    }

    Vector sumCols(Matrix mat){
        //returns the sum for each column.
        double outs [] = new double [mat.columns()];
        for (int i=0;i<mat.columns();i++){
            outs[i] = mat.getColumn(i).sum();
        }
        return Vector.fromArray(outs);
    }

    Matrix subtractCols(Matrix mat, Vector v){
        Matrix mat2 = mat.copy();
        //subtracts vector of values from each row. iterates by row
        for (int i=0;i<mat2.rows();i++){
            mat2.setRow(i,mat2.getRow(i).subtract(v));
        }
        return mat2;
    }

    Matrix integrateCols(Matrix mat){
        //Log.d("integrate mat size", String.valueOf(mat.rows())+ " " +String.valueOf(mat.columns()));
        Matrix out = mat.copy();
        Vector runSum = Vector.zero(out.columns());
        double fact = 2.00000;

        for (int i = 1; i< out.rows();i++){
            runSum = runSum.add(  mat.getRow(i).add(mat.getRow(i-1)).divide(fact)  );
            out.setRow(i, runSum);
            //Log.d("runsum",runSum.toString());
        }
        //Log.d("int out",out.toString());
        return out;
    }

    Matrix differentiateCols(Matrix mat){
        Matrix out = Matrix.zero(mat.rows(),mat.columns());
        for (int i = 1; i< mat.rows();i++){
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
        double toSubtract = 1.0000;
        double toDivide= 2.0000;
        Vector c = cross (A,B);
        Matrix added = Basic2DMatrix.identity(3).add(  ssc(c))  .  add( ssc(c).power(2).multiply  (  (toSubtract - A.hadamardProduct(B).sum()) / Math.pow(c.norm(),toDivide)  )  );
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