//much of source is adapted from https://github.com/ejoebstl/Android-Sensor-Log/blob/master/app/src/main/java/io/iam360/sensorlog/MainActivity.java

package capstone.wheelchairmonitor;


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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
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
                // initialize arrays
                ArrayList<Double> accX=new ArrayList<Double>();
                ArrayList<Double> accY=new ArrayList<Double>();
                ArrayList<Double> accZ=new ArrayList<Double>();
                ArrayList<Double> accT=new ArrayList<Double>();
                ArrayList<Double> magX=new ArrayList<Double>();
                ArrayList<Double> magY=new ArrayList<Double>();
                ArrayList<Double> magZ=new ArrayList<Double>();
                ArrayList<Double> magT=new ArrayList<Double>();
                //read in file
                if (toRead!=null) {
                    try {
                        Log.d("processing", "to read in not null");
                        BufferedReader br = new BufferedReader(new FileReader(toRead));
                        String line;


                        while ((line = br.readLine()) != null) {
                            //Log.d("processing line", line);
                            String[] values = line.split(";");
                            //extract vectors of measurements

                            if (values[1].equals("ACC")) {
                                accX.add(Double.parseDouble(values[2]));
                                accY.add(Double.parseDouble(values[3]));
                                accZ.add(Double.parseDouble(values[4]));
                                accT.add(Double.parseDouble(values[0])/(1.00e9)); //outputs in nanoseconds-divide by 10^9

                                continue;
                            }
                            if (values[1].equals("MAG")) {
                                magX.add(Double.parseDouble(values[2]));
                                magY.add(Double.parseDouble(values[3]));
                                magZ.add(Double.parseDouble(values[4]));
                                magT.add(Double.parseDouble(values[0])/(1.00e9)); //outputs in nanoseconds-divide by 10^9
                                continue;
                            }

                        }
                        br.close();
                    } catch (IOException e) {
                        //could not read file
                    }
                    //do processing.

                    //find effective sampling rates - acc is 500 hz and mag is 100 hz..
                    int nMag=magX.size();
                    Log.d("nMag", String.valueOf(nMag));
                    int nAcc=accX.size();
                    Log.d("nAcc", String.valueOf(nAcc));
                    Double len = accT.get(nAcc - 1) - accT.get(0);
                    Log.d("len", String.valueOf(len));
                    if  (nMag>0 & nAcc>0) {
                        double Fs_acc = nAcc / len;
                        Log.d("Fs_acc", String.valueOf(Fs_acc));
                        double Fs_mag = nMag / len;
                        Log.d("Fs_mag", String.valueOf(Fs_mag));
                    }
                    double sampratio = Math.floor(nAcc/nMag);
                    //downsample - done in naive way. very consistent 5 mag to 1 acc
                    for (int k = nAcc-1; k-- > 0; ) {
                        //Log.d("k", String.valueOf(k));
                        if ((k%sampratio)!=0) {
                            accX.remove(k);
                            accY.remove(k);
                            accZ.remove(k);
                            accT.remove(k);
                            continue;
                        }
                    }

                    //make lengths work together - no evidence for failure here yet!
                    while (nMag<accT.size()) {accX.remove(accT.size()-1);accY.remove(accT.size()-1);accZ.remove(accT.size()-1);accT.remove(accT.size()-1);}
                    while (magX.size()>accX.size()) {magX.remove(magT.size()-1);magY.remove(magT.size()-1);magZ.remove(magT.size()-1);magT.remove(magT.size()-1);}
                    nMag=magX.size();
                    Log.d("nMag down", String.valueOf(nMag));
                    nAcc=accX.size();
                    Log.d("nAcc down", String.valueOf(nAcc));

                    //decide which time vector to use, mag or acc. mag is less noisy (?). will make own idealized version
                    double Fs_effective = nAcc/len;
                    //tvec= //for loop

                    //remove gravity from acceleration

                    //do bandpass [0.3, 10] used in matlab

                    //integrate twice for position

                    //find frequency through peak finding

                    //calculate rotation area with polar approx

                    //TODO finish processing




                }



                return true;
            }
        });






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
                String[] values = line.split(";");
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

                String[] values = line.split(";");
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
                        writer.write(String.format("%d;ACC;%f;%f;%f;%f;%f;%f\n", evt.timestamp, evt.values[0], evt.values[1], evt.values[2], 0.f, 0.f, 0.f));
                        //Log.d("sensorwriting", "wrote acc\n");
                        break;
                    case Sensor.TYPE_MAGNETIC_FIELD:
                        if (!magFirst){magFirst = true;} //mag takes longer to spin up
                        //timestamp in nanoseconds
                        writer.write(String.format("%d;MAG;%f;%f;%f;%f;%f;%f\n", evt.timestamp, evt.values[0], evt.values[1], evt.values[2], 0.f, 0.f, 0.f));
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