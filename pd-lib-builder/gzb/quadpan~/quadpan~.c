/* code for "quadpan~.c" pd class.  This takes 1 signal and 4 arguments which determine the amount of the signal coming from 4 speakers.
 
 Speakers:
    1. Front Left
    2. Front Right
    3. Back Left
    4. Back Right
 
 Arguments:
    arg1 (left-right): 1 = full left, 0 = full right
    arg2 (right-left): reciprocal of arg 1 (change code compute!)
    arg3 (front-back): same logic
    arg4 (back-front): same logic
 
 To do: use only two argument one for left right (0 to 1) and one for front back (0 to 1)
 */

//github.com/pure-data/externals-howto#a-signal-external-pan
//maxpat source: [Somnitecs 4.1 panner v1.1] ->>> made my external from this patches logic.

#include "m_pd.h" //obligatory header
#include <math.h> //for the log e function

//must have
static t_class *quadpan_tilde_class;

typedef struct _quadpan_tilde {
    t_object x_obj; //must have
    
    //do it stupidly then optimize
    t_sample f_quadpan1; //mixing factor 1
    t_sample f_quadpan2; //mixing factor 2
    t_sample f_quadpan3; //mixing factor 3
    t_sample f_quadpan4; //mixing factor 4
    
    t_sample f;//believe this is for floats sent to signal in.
    
    //inlet one is for the signals to pan.
    t_inlet *x_in2;//inlet 2 (for messages)
    t_inlet *x_in3;//inlet 3 (for messages)
    t_inlet *x_in4;//inlet 4 (for messages)
    t_inlet *x_in5;//inlet 5 (for messages)
    
    t_outlet *x_out1;// outlet one (FL)
    t_outlet *x_out2;// outlet two (FR)
    t_outlet *x_out3;// outlet three (BL)
    t_outlet *x_out4;// outlet four (BR)
    
} t_quadpan_tilde; //must have

//perform routine:
//input arg: array with pointers
t_int *quadpan_tilde_perform(t_int *w)
{
    t_quadpan_tilde *x = (t_quadpan_tilde *)(w[1]);//object itself
    t_sample  *in1 =     (t_sample *)(w[2]);//inlet 1 (signal in1)
    t_sample  *out1 =    (t_sample *)(w[3]);//outlet (signal out1)
    t_sample  *out2 =    (t_sample *)(w[4]);//outlet (signal out2)
    t_sample  *out3 =    (t_sample *)(w[5]);//outlet (signal out3)
    t_sample  *out4 =    (t_sample *)(w[6]);//outlet (signal out4)
    int          n =            (int)(w[7]);//block size
    
    //limit the values of f_quadpanX between 0 and 1 (stupid way...can make it into function)
    t_sample f_quadpan1 =
    (x->f_quadpan1<0) ? 0.0 : (x->f_quadpan1>1) ? 1.0 : x->f_quadpan1;
    t_sample f_quadpan2 =
    (x->f_quadpan2<0) ? 0.0 : (x->f_quadpan2>1) ? 1.0 : x->f_quadpan2;
    t_sample f_quadpan3 =
    (x->f_quadpan3<0) ? 0.0 : (x->f_quadpan3>1) ? 1.0 : x->f_quadpan3;
    t_sample f_quadpan4 =
    (x->f_quadpan4<0) ? 0.0 : (x->f_quadpan4>1) ? 1.0 : x->f_quadpan4;

    //post("before log f_quadpan1 %f", f_quadpan1);
    
    //must do 4 times as values are different, the calculation could be made into a function though
    f_quadpan1 = (log(f_quadpan1*100 + 1))/4.62;
    f_quadpan2 = (log(f_quadpan2*100 + 1))/4.62;
    f_quadpan3 = (log(f_quadpan3*100 + 1))/4.62;
    f_quadpan4 = (log(f_quadpan4*100 + 1))/4.62;
    //post("after log f_quadpan1 %f", f_quadpan1);
    
    while (n--) {
        float samp = *in1++;//this fixed it somehow
        *out1++ = (samp)*(f_quadpan1*f_quadpan3);
        *out2++ = (samp)*(f_quadpan2*f_quadpan3);
        *out3++ = (samp)*(f_quadpan1*f_quadpan4);
        *out4++ = (samp)*(f_quadpan2*f_quadpan4);
    }

//    post("out1 gain: %f", f_quadpan1*f_quadpan3);
//    post("out2 gain: %f", f_quadpan2*f_quadpan3);
//    post("out3 gain: %f", f_quadpan1*f_quadpan4);
//    post("out4 gain: %f", f_quadpan2*f_quadpan4);

//    post("out1[n] %f", out1[n]);//post once every block
//    post("out2[n] %f", out2[n]);//post once every block
//    post("out3[n] %f", out3[n]);//post once every block
//    post("out4[n] %f", out4[n]);//post once every block

    //the return argument equals the argument of the perform-routine plus the number of pointer variables (as declared by the second argument of dsp_add) plus one
    return (w+8);
}

//signal processing method defined
//n = number of following pointer-arguments
void quadpan_tilde_dsp(t_quadpan_tilde *x, t_signal **sp)
{
    dsp_add(quadpan_tilde_perform, //following arguments passed to perform
            7,            //number of pointers
            x,            //object itself
            sp[0]->s_vec, //sp[0] is the first in-signal,
            sp[1]->s_vec, //sp[1] points to the out-signal 1 (FL)
            sp[2]->s_vec, //sp[2] points to the out-signal 2 (FR)
            sp[3]->s_vec, //sp[3] points to the out-signal 3 (BL)
            sp[4]->s_vec, //sp[4] points to the out-signal 4 (BR)
            sp[0]->s_n);  //length of first in-signal (block size = N)
}

//destrutor
void quadpan_tilde_free(t_quadpan_tilde *x)
{
    //clear inlets
    inlet_free(x->x_in2);
    inlet_free(x->x_in3);
    inlet_free(x->x_in4);
    inlet_free(x->x_in5);
    
    //clear outlets
    outlet_free(x->x_out1);
    outlet_free(x->x_out2);
    outlet_free(x->x_out3);
    outlet_free(x->x_out4);
}

//creator
void *quadpan_tilde_new(t_floatarg f, t_floatarg g, t_floatarg h, t_floatarg i) {
    
    //pass class from setup to default method and cast
    t_quadpan_tilde *x = (t_quadpan_tilde *)pd_new(quadpan_tilde_class);
    
    //store value in obj data space
    x->f_quadpan1 = f; //post("quadpan1 inlet %f", f)
    x->f_quadpan2 = g;
    x->f_quadpan3 = h;
    x->f_quadpan4 = i;
    
    //make float inlets, inlet one reserved for signal in
    x->x_in2 = floatinlet_new (&x->x_obj, &x->f_quadpan1);
    x->x_in3 = floatinlet_new (&x->x_obj, &x->f_quadpan2);
    x->x_in4 = floatinlet_new (&x->x_obj, &x->f_quadpan3);
    x->x_in5 = floatinlet_new (&x->x_obj, &x->f_quadpan4);
    
    //make signal outlets (FL, FR, BL, BR)
    x->x_out1 = outlet_new(&x->x_obj, &s_signal);
    x->x_out2 = outlet_new(&x->x_obj, &s_signal);
    x->x_out3 = outlet_new(&x->x_obj, &s_signal);
    x->x_out4 = outlet_new(&x->x_obj, &s_signal);
    
    //return casted pointer
    return (void *)x;
    
}

//SETUP
void quadpan_tilde_setup(void) {
    
    quadpan_tilde_class =   class_new(gensym("quadpan~"),//standard
                                  (t_newmethod)quadpan_tilde_new,//creator
                                  (t_method)quadpan_tilde_free,//destructor
                                  sizeof(t_quadpan_tilde),//allocation
                                  CLASS_DEFAULT,//appearance
                                  A_DEFFLOAT,//creation arg (LR) ->1=LEFT
                                  A_DEFFLOAT,//creation arg (RL) ->1=RIGHT
                                  A_DEFFLOAT,//creation arg (FB) ->1=FRONT
                                  A_DEFFLOAT,//creation arg (BF) ->1=BACK
                                  0);//termination
    
    //called when DSP is turned on
    class_addmethod(quadpan_tilde_class,//specify class to add method to
                    (t_method)quadpan_tilde_dsp,//cast method name to pd type
                    gensym("dsp"),//id for dsp objects
                    A_CANT,//prevent crashes
                    0);//termination
    
    //makes inlet one a signal inlet
    CLASS_MAINSIGNALIN(quadpan_tilde_class, t_quadpan_tilde, f);
    //third argument of CLASS_MAINSIGNALIN is a (dummy-)floating point-variable of the data space, that is needed to automatically convert “float”-messages into signals if no signal is present at the signal-inlet.
}

