/* code for "foaAmbiPan~.c" pd class.  This takes 1 signal and 2 arguments which determine the elevation and azimuth. It also has 4 ouputs since we need to get W, X, Y and Z.
 
 Arguments:
    arg1 (azimuth): 180 to -180 where 0 is center, 90 is left and -90 is right.
    arg2 (elevation): 180 to -180 where 0 is flat, 90 is directly above and -90 is directly below
 
 Eq: to simplify our code k will always be equal to 1
    W: 1/k * s (1/sqrt(2))
    X: 1/k * s(cos(a)cos(e))
    Y: 1/k * s(sin(a)cos(e))
    Z: 1/k * s(sin(e))
 
 where k is the number of signals, a is azimuth and e is elevation.
 
 Notes: we have to decide which ordering we will use. I think WXYZ is less confusing, that one is called Furse-Malham. WYZX is ACN. We have to be careful about radians v. degrees!
 */

#include "m_pd.h" //obligatory header
#include <math.h> //for the math stuff (sin, cos, ...)

//must have
static t_class *foaAmbiPan_tilde_class;

//this is essentially the architecture of the object, aka a sort of class definition.
typedef struct _foaAmbiPan_tilde {
    
    t_object x_obj; //must have, it's like a reference to itself
    
    //do it stupidly then optimize
    t_sample f_azimuth; //left to right
    t_sample f_elevation; //up to down
    
    t_sample f;//believe this is for floats sent to signal in.
    
    //inlet one is for the signal that we are going to pan
    t_inlet *x_in2;//inlet 2 (for messages)
    t_inlet *x_in3;//inlet 3 (for messages)
    
    t_outlet *x_out1;// outlet one (W)
    t_outlet *x_out2;// outlet two (X)
    t_outlet *x_out3;// outlet three (Y)
    t_outlet *x_out4;// outlet four (Z)
    
} t_foaAmbiPan_tilde; //must have

//perform routine:
//input arg: array with pointers
t_int *foaAmbiPan_tilde_perform(t_int *w)
{
    t_foaAmbiPan_tilde *x = (t_foaAmbiPan_tilde *)(w[1]);//object itself
    t_sample  *in1  =    (t_sample *)(w[2]);//inlet 1 (signal in1)
    t_sample  *out1 =    (t_sample *)(w[3]);//outlet (signal out1)
    t_sample  *out2 =    (t_sample *)(w[4]);//outlet (signal out2)
    t_sample  *out3 =    (t_sample *)(w[5]);//outlet (signal out3)
    t_sample  *out4 =    (t_sample *)(w[6]);//outlet (signal out4)
    int          n =            (int)(w[7]);//block size
    
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
    //W: 1/k * s (1/sqrt(2))
    //X: 1/k * s(cos(a)cos(e))
    //Y: 1/k * s(sin(a)cos(e))
    //Z: 1/k * s(sin(e))
    
    while (n--) {
        //grab a single sample from the input
        float samp = *in1++;//this fixed it somehow
        
        *out1++ = (samp)*(1/(sqrt(2)));             //w
        *out2++ = (samp)*(cos(a)*cos(e));            //x
        *out3++ = (samp)*(sin(a)*cos(e));            //y
        *out4++ = (samp)*(sin(e));                  //z
    }

    //the return argument equals the argument of the perform-routine plus the number of pointer variables (as declared by the second argument of dsp_add) plus one
    return (w+8);
}

//signal processing method defined
//n = number of following pointer-arguments
void foaAmbiPan_tilde_dsp(t_foaAmbiPan_tilde *x, t_signal **sp)
{
    dsp_add(foaAmbiPan_tilde_perform, //following arguments passed to perform
            7,            //number of pointers
            x,            //object itself
            sp[0]->s_vec, //sp[0] is the first in-signal,
            sp[1]->s_vec, //sp[1] points to the out-signal 1 (W)
            sp[2]->s_vec, //sp[2] points to the out-signal 2 (X)
            sp[3]->s_vec, //sp[3] points to the out-signal 3 (Y)
            sp[4]->s_vec, //sp[4] points to the out-signal 4 (Z)
            sp[0]->s_n);  //length of first in-signal (block size = N)
}

//destrutor
void foaAmbiPan_tilde_free(t_foaAmbiPan_tilde *x)
{
    //clear inlets
    inlet_free(x->x_in2); //used for azimuth
    inlet_free(x->x_in3); //used for elevation
    
    //clear outlets
    outlet_free(x->x_out1); //output W
    outlet_free(x->x_out2); //output X
    outlet_free(x->x_out3); //output Y
    outlet_free(x->x_out4); //output Z
}

//creator
void *foaAmbiPan_tilde_new(t_floatarg f, t_floatarg g) {
    
    //pass class from setup to default method and cast
    t_foaAmbiPan_tilde *x = (t_foaAmbiPan_tilde *)pd_new(foaAmbiPan_tilde_class);
    
    //store value in obj data space
    //post("quadpan1 inlet %f", f)
    x->f_azimuth = f;
    x->f_elevation = g;
    
    //make float inlets, inlet one reserved for signal in
    //args(data space, specific address to write values to)
    x->x_in2 = floatinlet_new (&x->x_obj, &x->f_azimuth);
    x->x_in3 = floatinlet_new (&x->x_obj, &x->f_elevation);
    
    //make signal outlets (W, X, Y, Z)
    x->x_out1 = outlet_new(&x->x_obj, &s_signal);
    x->x_out2 = outlet_new(&x->x_obj, &s_signal);
    x->x_out3 = outlet_new(&x->x_obj, &s_signal);
    x->x_out4 = outlet_new(&x->x_obj, &s_signal);
    
    //return casted pointer
    return (void *)x;
}

//SETUP
void foaAmbiPan_tilde_setup(void) {
    
    foaAmbiPan_tilde_class =   class_new(gensym("foaAmbiPan~"),//standard
                                  (t_newmethod)foaAmbiPan_tilde_new,//creator
                                  (t_method)foaAmbiPan_tilde_free,//destructor
                                  sizeof(t_foaAmbiPan_tilde),//allocation
                                  CLASS_DEFAULT,//appearance
                                  A_DEFFLOAT,//creation arg (azimuth)
                                  A_DEFFLOAT,//creation arg (elevation)
                                  0);//termination
    
    //called when DSP is turned on
    class_addmethod(foaAmbiPan_tilde_class,//specify class to add method to
                    (t_method)foaAmbiPan_tilde_dsp,//cast method name to pd type
                    gensym("dsp"),//id for dsp objects
                    A_CANT,//prevent crashes
                    0);//termination
    
    //makes inlet one a signal inlet
    CLASS_MAINSIGNALIN(foaAmbiPan_tilde_class, t_foaAmbiPan_tilde, f);
    //third argument of CLASS_MAINSIGNALIN is a (dummy-)floating point-variable of the data space, that is needed to automatically convert “float”-messages into signals if no signal is present at the signal-inlet.
}

