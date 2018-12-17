/* code for "foaRot~.c" pd class. This external takes 2 arguments which determine the elevation and azimuth and 4 signals (W, X, Y, Z) (B-format FOA signals). It also has 4 ouputs since we need to get W, X, Y and Z.

 Arguments:
    arg1 (azimuth): 180 to -180 where 0 is center, 90 is left and -90 is right.
    arg2 (elevation): 180 to -180 where 0 is flat, 90 is directly above and -90 is directly below

 Eq:    yaw (Z) : X' = X * cos(a) - Y * sin(a) & Y' = X * sin(a) + Y * cos(a)
        pitch(Y) : X' = X * cos(b) - Z * sin(b) & Z' = X * sin(b) + Z * cos(b)
        roll (X) : Y’ = Y * cos(b) – Z * sin(b) & Z’ = Y * sin(b) + Z * cos(b)

 where a is the azimuth angle and b is the elevation angle. Head trackers send our three values however (X,Y,Z), I am unclear as to how those values fit this model.

 Notes: we have to decide which ordering we will use. I think WXYZ is less confusing, that one is called Furse-Malham. WYZX is ACN. We have to be careful about radians v. degrees! TODO (multiply by pi/180)
 */

#include "m_pd.h" //obligatory header
#include <math.h> //for the math stuff (sin, cos, ...)

//must have
static t_class *foaRot_tilde_class;

//this is essentially the architecture of the object, aka a sort of class definition.
typedef struct _foaRot_tilde {

    t_object x_obj; //must have, it's like a reference to itself

    //do it stupidly then optimize
    t_sample f_azimuth; //left to right
    t_sample f_elevation; //up to down

    t_sample f;//believe this is for floats sent to main signal in.

    //we need 4 signal inlets and 2 float inlets

    //      *x_in1 = MAINSIGNALIN
    t_inlet *x_in2; // X
    t_inlet *x_in3; // Y
    t_inlet *x_in4; // Z

    //azimuth and elevation inlets
    t_inlet *x_in5;//inlet 2 (for messages) [azimuth]
    t_inlet *x_in6;//inlet 3 (for messages) [elevation]

    t_outlet *x_out1;// outlet one (W)
    t_outlet *x_out2;// outlet two (X)
    t_outlet *x_out3;// outlet three (Y)
    t_outlet *x_out4;// outlet four (Z)

} t_foaRot_tilde; //must have

//////////////////
//PERFORM ROUTINE:
//////////////////

//input arg: array with pointers
t_int *foaRot_tilde_perform(t_int *w)
{
    //i dont remember what 0 is but i think its a pointer to a pointer to the object?
    t_foaRot_tilde *x =  (t_foaRot_tilde *)(w[1]);//object itself
    t_sample  *in1  =    (t_sample *)(w[2]);//inlet 1 (signal in1)
    t_sample  *in2  =    (t_sample *)(w[3]);//inlet 1 (signal in1)
    t_sample  *in3  =    (t_sample *)(w[4]);//inlet 1 (signal in1)
    t_sample  *in4  =    (t_sample *)(w[5]);//inlet 1 (signal in1)
    t_sample  *out1 =    (t_sample *)(w[6]);//outlet (signal out1)
    t_sample  *out2 =    (t_sample *)(w[7]);//outlet (signal out2)
    t_sample  *out3 =    (t_sample *)(w[8]);//outlet (signal out3)
    t_sample  *out4 =    (t_sample *)(w[9]);//outlet (signal out4)
    int          n =            (int)(w[10]);//block size

    //limit the values of elevation and azimuth between -180 and 180
    //TODO
    //        t_sample azimuth   = x->azimuth % 180;
    //        t_sample elevation = x->elevation % 180;

    //convert degrees to radians by multiplying by pi/180. M_PI defined in math.h
    t_sample a = x->f_azimuth * (M_PI/180);
    t_sample e = x->f_elevation * (M_PI/180);

//    post("azimuth = %f", a);
//    post("elevation = %f", e);

    //the equation, for reference
//    yaw (Zrot)  : X' = X * cos(a) - Y * sin(a) & Y' = X * sin(a) + Y * cos(a)
//    pitch(Yrot) : X' = X * cos(b) - Z * sin(b) & Z' = X * sin(b) + Z * cos(b)
//    roll (Xrot) : Y’ = Y * cos(b) – Z * sin(b) & Z’ = Y * sin(b) + Z * cos(b)

    while (n--) {
        //grab a single sample from each input
        float w_samp = *in1++;//this fixed it somehow
        float x_samp = *in2++;//this fixed it somehow
        float y_samp = *in3++;//this fixed it somehow
        float z_samp = *in4++;//this fixed it somehow

        //W always stays the same

        //X Out
        *out2++ =   ( (x_samp * cos(a)) - (y_samp * sin(a)) + (x_samp * cos(e)) - (z_samp * sin(e)));

        //Y Out
        *out3++ =   ( (y_samp * cos(e)) - (z_samp * sin(e)) + (x_samp * sin(a)) + (y_samp * cos(a)));

        //Z Out
        *out4++ =   ( (x_samp * sin(e)) + (z_samp * cos(e)) + (y_samp * sin(e)) + (z_samp * cos(e)));
    }

    //the return argument equals the argument of the perform-routine plus the number of pointer variables (as declared by the second argument of dsp_add) plus one
    return (w+8);
}

//signal processing method defined
//n = number of following pointer-arguments
void foaRot_tilde_dsp(t_foaRot_tilde *x, t_signal **sp)
{
    dsp_add(foaRot_tilde_perform, //following arguments passed to perform
            10,            //number of pointers [ten total?]
            x,            //object itself
            sp[0]->s_vec, //sp[0] is the first in-signal,   [Wi]
            sp[1]->s_vec, //sp[1] is the 2nd   in-signal    [Xi]
            sp[2]->s_vec, //sp[2] is the 3rd   in-signal    [Yi]
            sp[3]->s_vec, //sp[3] is the 4th   in-signal    [Zi]
            sp[4]->s_vec, //sp[4] is the first out-signal   [Wo]
            sp[5]->s_vec, //sp[5] is the 2nd   out-signal   [Xo]
            sp[6]->s_vec, //sp[5] is the 3rd   out-signal   [Yo]
            sp[7]->s_vec, //sp[7] is the 4th   out-signal   [Zo]
            sp[0]->s_n);  //length of first in-signal (block size = N)
}

//destrutor
void foaRot_tilde_free(t_foaRot_tilde *x)
{
    //clear inlets
    inlet_free(x->x_in2); //used for X
    inlet_free(x->x_in3); //used for Y
    inlet_free(x->x_in4); //used for Z

    //clear outlets
    outlet_free(x->x_out1); //output W
    outlet_free(x->x_out2); //output X
    outlet_free(x->x_out3); //output Y
    outlet_free(x->x_out4); //output Z
}

//creator
void *foaRot_tilde_new(t_floatarg f, t_floatarg g) {

    //pass class from setup to default method and cast
    t_foaRot_tilde *x = (t_foaRot_tilde *)pd_new(foaRot_tilde_class);

    //store value in obj data space
    //post("quadpan1 inlet %f", f)
    x->f_azimuth = f;
    x->f_elevation = g;

    //initialize signal inlets
    //MAINSIGNALIN does not need to be initialized
    x->x_in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);//X
    x->x_in3 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);//Y
    x->x_in4 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);//Z

    //make float inlets
    //args(data space, specific address to write values to)
    x->x_in5 = floatinlet_new (&x->x_obj, &x->f_azimuth);
    x->x_in6 = floatinlet_new (&x->x_obj, &x->f_elevation);

    //make signal outlets (W, X, Y, Z)
    x->x_out1 = outlet_new(&x->x_obj, &s_signal);
    x->x_out2 = outlet_new(&x->x_obj, &s_signal);
    x->x_out3 = outlet_new(&x->x_obj, &s_signal);
    x->x_out4 = outlet_new(&x->x_obj, &s_signal);

    //return casted pointer
    return (void *)x;
}

//SETUP
void foaRot_tilde_setup(void) {

    foaRot_tilde_class =   class_new(gensym("foaRot~"),//standard
                                  (t_newmethod)foaRot_tilde_new,//creator
                                  (t_method)foaRot_tilde_free,//destructor
                                  sizeof(t_foaRot_tilde),//allocation
                                  CLASS_DEFAULT,//appearance
                                  A_DEFFLOAT,//creation arg (azimuth)
                                  A_DEFFLOAT,//creation arg (elevation)
                                  0);//termination

    //called when DSP is turned on
    class_addmethod(foaRot_tilde_class,//specify class to add method to
                    (t_method)foaRot_tilde_dsp,//cast method name to pd type
                    gensym("dsp"),//id for dsp objects
                    A_CANT,//prevent crashes
                    0);//termination

    //makes inlet one a signal inlet with extra functionality...
    CLASS_MAINSIGNALIN(foaRot_tilde_class, t_foaRot_tilde, f);
    //third argument of CLASS_MAINSIGNALIN is a (dummy-)floating point-variable of the data space, that is needed to automatically convert “float”-messages into signals if no signal is present at the signal-inlet.
}
