#include "m_pd.h"
#include <stdlib.h>
#include <math.h>

#ifdef NT //idk what this is, TRE code.
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

/* ------------------------ pool_vox~ ----------------------------- */

/* Code for pool_vox~. This one voice out of 16 which will be used to make the new version of Pool. Eventually I want to get 64 voices but for now we will do 16.
 
 Arguments:
 arg1: hue, used to determine frequency of oscillator [square wave]
 arg2: saturation, used to determine fc of lp filter (one-pole)
 arg3: value, used to determine volume of oscillator
 the value will be 0 when the image is black, we will do 1 minus value instead
    this will give us sound only when there is color (white = no sound)
 arg4: azimuth, will determine the position of the voice in 3OA (horizontal)
 arg5: elevation, will determine the position of the voice in 3OA (vertical)
 
 outlets: 16 ambisonic channels of encoded audio. SN3D (ACN)
 */

#define WAVETABLESIZE 1024

static t_class *pool_vox_class;

typedef struct _pool_vox
{
    t_object x_obj; 	/* obligatory header */
    t_float x_f;    	/* place to hold inlet's value if it's set by message */
    t_float *wavetable;	/* a place to hold the square wave */
    t_float phase;
    t_float samplerate;
    t_float hue;        // hue from GEM
    //  (x_f = hue)
    t_float sat;        // saturation from GEM
    t_float value;      // value or luminance from GEM
    t_float azi;        // azimuth of this voice
    t_float elev;       // elevation of this voice
    
    //stuff for filter, controlled by saturation from GEM
    /*---------------------------------------------------*/
    t_float cutoff; //filter cutoff modulated by sat
    t_float filt_q; //filter quality or q
    t_float v0, v1, temp; //hold sample values for filter.
    /*---------------------------------------------------*/
    
    //inlet one is for the signal (hue: frequency of square wave)
    t_inlet *x_in2;//inlet 2 (for messages - saturation: fc)
    t_inlet *x_in3;//inlet 3 (for messages - value: volume of voice)
    
    //decided to use creation arg instead
    //t_inlet *x_in4;//inlet 4 (for messages - azimuth)
    //t_inlet *x_in5;//inlet 5 (for messages - elevation)
    
    //16 outlets for 3OA
    //note there will be 16 voices, each with 16 ambi channels. 256 total signals.
    //heavy computation.
    t_outlet *x_out1, *x_out2, *x_out3, *x_out4, *x_out5, *x_out6, *x_out7, *x_out8, *x_out9, *x_out10, *x_out11, *x_out12, *x_out13, *x_out14, *x_out15, *x_out16;
    
    //coefficients for amplitude based on 3OA (16 channels)
    t_float W, Y, Z, X, V, T, R, S, U, Q, O, M, K, L, N, P;
    
} t_pool_vox;

// some nice little interpolation routines
// no interp [1]
static inline float no_interpolate(t_pool_vox *x)
{
    int x_1 = x->phase * WAVETABLESIZE;
    return(x->wavetable[x_1 % WAVETABLESIZE]);
}
//linear interpolation [2]
static inline float lin_interpolate(t_pool_vox *x)
{
    int x_1 = x->phase * WAVETABLESIZE;
    float y_1 = x->wavetable[x_1 % WAVETABLESIZE];
    float y_2 = x->wavetable[(x_1 + 1) % WAVETABLESIZE];
    return (y_2 - y_1) * ((x->phase * WAVETABLESIZE) - x_1) + y_1;
}
//using quad interp (hermitian) [3]
static inline float quad_interpolate(t_pool_vox *x)
{
    int truncphase = (int) (x->phase * WAVETABLESIZE);
    float fr = (x->phase * WAVETABLESIZE) - ((float) truncphase);
    float inm1 = x->wavetable[(truncphase - 1) % WAVETABLESIZE];
    float in   = x->wavetable[(truncphase + 0) % WAVETABLESIZE];
    float inp1 = x->wavetable[(truncphase + 1) % WAVETABLESIZE];
    float inp2 = x->wavetable[(truncphase + 2) % WAVETABLESIZE];
    
    return in + 0.5 * fr *
    (inp1 - inm1 + fr *
        (4.0 * inp1 + 2.0 * inm1 - 5.0 * in - inp2 + fr *
            (3.0 *
                (in - inp1)
                    - inm1 + inp2)
                                    )
                                        );
}

//static inline float two_pole_lp(t_pool_vox *x)
//{
//
//    int truncphase = (int) (x->phase * WAVETABLESIZE);
//    float fr = (x->phase * WAVETABLESIZE) - ((float) truncphase);
////    float inm1 = x->wavetable[(truncphase - 1) % WAVETABLESIZE];
//    float in   = x->wavetable[(truncphase + 0) % WAVETABLESIZE];
////    float inp1 = x->wavetable[(truncphase + 1) % WAVETABLESIZE];
////    float inp2 = x->wavetable[(truncphase + 2) % WAVETABLESIZE];
//    //return = v1;
//}

/* this is the actual performance routine which acts on the samples.
 It's called with a single pointer "w" which is our location in the
 DSP call list.  We return a new "w" which will point to the next item
 after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
 (no use to us), w[1] and w[2] are the input and output vector locations,
 and w[3] is the number of points to calculate. */

static t_int *pool_vox_perform(t_int *w)
{
    //dsp routine, only deals with signal ins/outs
    //  t_floats are temporarily used to store sample values for outputs
    t_pool_vox *x = (t_pool_vox *)(w[1]); //self reference
    t_float *freq = (t_float *)(w[2]); //first inlet (hue)  [float pointer]
    
    //sixteen ambisonic channels
    t_float *out1 = (t_float *)(w[3]); t_float *out2 = (t_float *)(w[4]);
    t_float *out3 = (t_float *)(w[5]); t_float *out4 = (t_float *)(w[6]);
    t_float *out5 = (t_float *)(w[7]); t_float *out6 = (t_float *)(w[8]);
    t_float *out7 = (t_float *)(w[9]); t_float *out8 = (t_float *)(w[10]);
    t_float *out9 = (t_float *)(w[11]); t_float *out10 = (t_float *)(w[12]);
    t_float *out11 = (t_float *)(w[13]); t_float *out12 = (t_float *)(w[14]);
    t_float *out13 = (t_float *)(w[15]); t_float *out14 = (t_float *)(w[16]);
    t_float *out15 = (t_float *)(w[17]); t_float *out16 = (t_float *)(w[18]);
    
    int n = (int)(w[19]); //block size
    
    //get saturation ad value from main struct
    t_sample s_sat = x->sat;        //used for cutoff of 2-pole filt.
    t_sample s_value = x->value;    //1-value = amplitude of voice
    
    // i like counting from zero
    int blocksize = n;
    int i, sample = 0;      //sample here is used to index, not store amplitudes
    float phaseincrement;
    float findex;
    float wphase;
    int	iindex;
    float samp;             //samp is the actual value
    float f_res;            //resonance of 2-pole filter
    float f_cutoff;         //cutoff storage
    
    //resonance fixed to 0.0625, filt_q set to 16
    f_res = 1.0f/x->filt_q;
    
    while (n--)
    {
        // first we need to calculate the phase increment from the frequency
        // and sample rate - this is the number of cycles per sample
        // freq = cyc/sec, sr = samp/sec, phaseinc = cyc/samp = freq/sr
        phaseincrement = *(freq+sample)/x->samplerate;
        
        // now, increment the phase and make sure it doesn't go over 1.0
        x->phase += phaseincrement;
        while(x->phase >= 1.0f)
            x->phase -= 1.0f;
        while(x->phase < 0.0f)
            x->phase += 1.0f;
        wphase = x->phase * x->phase;
        // to test - lets just spit out phase
        *(out1+sample) = x->phase;
        // now grab the sample from the table
        findex = WAVETABLESIZE * wphase;
        iindex = (int)findex;
        //		*(out1+sample) = *(x->wavetable + iindex);
        //		*(out1+sample) = lin_interpolate(x);
        
        /*
         W                  0th order (0)
         Y Z X              1st order (1, 2, 3)
         V T R S U          2nd order (4, 5, 6, 7, 8)
         Q O M K L N P      3rd order (9, 10, 11, 12, 13, 14, 15)
         */
        
        //after you get the sample from the wavetable, interpolate.
        //  it wasn't working, just had to restart Pd to update external. :)
        samp = quad_interpolate(x);
        
        //multiply samp by 1 minus "value" from GEM [HSV].
        samp = samp * (1.0f - s_value);

        //filter w/ cutoff based on saturation from GEM.
        // range 0-1 from saturation multiplied by 256
        
        //different tests. checking out what sounds good.
        //x->cutoff = 2.0f * sinf((float)M_PI * (s_sat * (x->samplerate / 4.0f)) / x->samplerate);
        //x->cutoff = 2.0f * sinf((float)M_PI * (s_sat * (x->samplerate/16.0f)) / x->samplerate);
        
        //256 = 16 x 16, arbitrary [settled on this]
        //to-do: nicer sounding filter
        x->cutoff = 2.0f * sinf((float)M_PI * (s_sat * (256.0f)) / x->samplerate);
        
        //mtof(64), 16X4
        //x->cutoff = 2.0f * sinf((float)M_PI * (s_sat * (329.628f)) / x->samplerate);

        
        //www.musicdsp.org/en/latest/Filters/185-1-rc-and-c-filter.html
        //implementing the filter
        x->v0 = (1.0f - f_res * x->cutoff) * x->v0 -
                x->cutoff * x->v1 +
                x->cutoff * samp;
        
        x->v1 = (1.0f - f_res * x->cutoff) * x->v1 +
                x->cutoff * x->v0;
        
        samp = x->v1;
        //end of filter implementation (to-do: x->temp not used)
        
        //then apply coeffs based on SN3D calculations
        *(out1+sample) = samp * x->W;
        *(out2+sample) = samp * x->Y;
        *(out3+sample) = samp * x->Z;
        *(out4+sample) = samp * x->X;
        *(out5+sample) = samp * x->V;
        *(out6+sample) = samp * x->T;
        *(out7+sample) = samp * x->R;
        *(out8+sample) = samp * x->S;
        *(out9+sample) = samp * x->U;
        *(out10+sample) = samp * x->Q;
        *(out11+sample) = samp * x->O;
        *(out12+sample) = samp * x->M;
        *(out13+sample) = samp * x->K;
        *(out14+sample) = samp * x->L;
        *(out15+sample) = samp * x->N;
        *(out16+sample) = samp * x->P;
        
        sample++;
    }
    //the return argument equals the argument of the perform-routine plus the number of pointer variables (as specified by the second argument of the dsp_add method) plus one. used to shift pointers [binary].
    return (w+20);
}

/* called to start DSP.  Here we call Pd back to add our perform
 routine to a linear callback list which Pd in turn calls to grind
 out the samples. */
static void pool_vox_dsp(t_pool_vox *x, t_signal **sp)
{
    // we'll initialize samplerate when starting up
    x->samplerate = sp[0]->s_sr;
    
    //passing w to pool_vox_perform
    dsp_add(pool_vox_perform,   //following args passed to perform
            19,                 //num pointers
            x,                  //object itself [1]
            sp[0]->s_vec,       //first input signal (hue/freq) [2]
            sp[1]->s_vec,       //out signal W [3]
            sp[2]->s_vec,       //out signal Y [4]
            sp[3]->s_vec,       //out signal Z [5]
            sp[4]->s_vec,       //out signal X [6]
            sp[5]->s_vec,       //out signal V [7]
            sp[6]->s_vec,       //out signal T [8]
            sp[7]->s_vec,       //out signal R [9]
            sp[8]->s_vec,       //out signal S [10]
            sp[9]->s_vec,       //out signal U [11]
            sp[10]->s_vec,      //out signal Q [12]
            sp[11]->s_vec,      //out signal O [13]
            sp[12]->s_vec,      //out signal M [14]
            sp[13]->s_vec,      //out signal K [15]
            sp[14]->s_vec,      //out signal L [16]
            sp[15]->s_vec,      //out signal N [17]
            sp[16]->s_vec,      //out signal P [18]
            sp[0]->s_n);        //len of 1st in signal, blocksize = N [19]
}

//constructor method
static void *pool_vox_new(t_floatarg f, t_floatarg g)
{
    //a for azi, e for elev.
    float twopi, size, a, e;
    int i, j;
    
    //allocate memory for the object itself
    t_pool_vox *x = (t_pool_vox *)pd_new(pool_vox_class);
    
    //pass creation arguments to struct
    //convert degrees to radians by multiplying by pi/180. M_PI defined in math.h
    a = x->azi = f * ((float)M_PI/180.0f);
    e = x->elev = g * ((float)M_PI/180.0f);
    
    //calculate 3OA coefficients (once during creation)
    //  source: blue ripple sound (SN3D, so ACN ordering)
    //  using angle/elevation representation,
    //      sqrt included in math.h
    /*-----------------------------------------------------------------------------*/
    //ACN Channel Number
    x->W = 1.0f;                                                //0
    x->Y = sin(a)*cos(e);                                       //1
    x->Z = sin(e);                                              //2
    x->X = cos(a)*cos(e);                                       //3
    x->V = sqrt(0.75f)*sin(2.0f*a)*cos(e)*cos(e);               //4
    x->T = sqrt(0.75f)*sin(a)*sin(2.0f*e);                      //5
    x->R = (0.5f)*(3.0f*sin(e)*sin(e)-1.0f);                    //6
    x->S = sqrt(0.75f)*cos(a)*sin(2.0f*e);                      //7
    x->U = sqrt(0.75f)*cos(2.0f*a)*cos(e)*cos(e);               //8
    x->Q = sqrt(0.625f)*sin(3.0f*a)*cos(e)*cos(e)*cos(e);       //9
    x->O = sqrt(3.75f)*sin(2.0f*a)*sin(e)*cos(e)*cos(e);        //10
    x->M = sqrt(0.375f)*sin(a)*cos(e)*(5.0f*sin(e)*sin(e)-1.0f);//11
    x->K = (0.5f)*sin(e)*(5.0f*sin(e)*sin(e)-3.0f);             //12
    x->L = sqrt(0.375f)*cos(a)*cos(e)*(5.0f*sin(e)*sin(e)-1.0f);//13
    x->N = sqrt(3.75f)*cos(2.0f*a)*sin(e)*cos(e)*cos(e);        //14
    x->P = sqrt(0.625f)*cos(3.0f*a)*cos(e)*cos(e)*cos(e);       //15
    /*-----------------------------------------------------------------------------*/
    
    //make float inlets, inlet 1 reserved for signal in (freq/hue)
    //  args(dataspace, specific address to write values to)
    //  note: azi and elev can't be changed once set, this is intentional.
    x->x_in2 = floatinlet_new (&x->x_obj, &x->sat);
    x->x_in3 = floatinlet_new (&x->x_obj, &x->value);
    
    //make signal outlets (16 total)
    x->x_out1 = outlet_new(&x->x_obj, &s_signal); x->x_out2 = outlet_new(&x->x_obj, &s_signal);
    x->x_out3 = outlet_new(&x->x_obj, &s_signal); x->x_out4 = outlet_new(&x->x_obj, &s_signal);
    x->x_out5 = outlet_new(&x->x_obj, &s_signal); x->x_out6 = outlet_new(&x->x_obj, &s_signal);
    x->x_out7 = outlet_new(&x->x_obj, &s_signal); x->x_out8 = outlet_new(&x->x_obj, &s_signal);
    x->x_out9 = outlet_new(&x->x_obj, &s_signal); x->x_out10 = outlet_new(&x->x_obj, &s_signal);
    x->x_out11 = outlet_new(&x->x_obj, &s_signal); x->x_out12 = outlet_new(&x->x_obj, &s_signal);
    x->x_out13 = outlet_new(&x->x_obj, &s_signal); x->x_out14 = outlet_new(&x->x_obj, &s_signal);
    x->x_out15 = outlet_new(&x->x_obj, &s_signal); x->x_out16 = outlet_new(&x->x_obj, &s_signal);
    
    // initialize variables to 0.0f [except azi and elev.]
    x->x_f = 0.0f;
    x->phase = 0.0f;
    
    x->sat = 0.0f;      //used to control cutoff of 2-pole LP filt.
    x->value = 0.0f;    //used to control amplitude of voice
    
    /*--------------2-pole filter stuff-----------------------------------------*/
    x->cutoff = 0.0f;   //2*sin(pi*f/fs) where f is controlled by saturation
        // saturation range 0-1 gets multiplied by 256 (arbitrary)
    x->filt_q = 16.0f;   //setting this to something fixed.
        //r = 1/filt_q (resonance of LP filter)
    x->v0 = 0.0f;       //storage for 2-pole filter (value 0)
    x->v1 = 0.0f;       //storage for 2-pole filter (value 1)
    x->temp = 0.0f;     //storage for 2-pole filter
    /*--------------2-pole filter stuff-----------------------------------------*/

    twopi = 8.0f * atanf(1.0f);
    size = (float)WAVETABLESIZE;
    
    // space for WAVETABLESIZE samples
    x->wavetable = (t_float *)malloc(WAVETABLESIZE * sizeof(t_float));
    
    // fill it up with a sine wave
    for(i = 0; i < WAVETABLESIZE; i++)
    {
        //16 harmonics.
//        for(j = 0; j < 15; j++)
//            {
            //*(x->wavetable+i) = 1.0f/(1 + j*2.0f) * sinf(twopi * (float)i/size);
                *(x->wavetable+i)  = 1.0f/3.0f * sinf(twopi * 1.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/3.0f * sinf(twopi * 3.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/5.0f * sinf(twopi * 5.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/7.0f * sinf(twopi * 7.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/9.0f * sinf(twopi * 9.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/11.0f * sinf(twopi * 11.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/13.0f * sinf(twopi * 13.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/15.0f * sinf(twopi * 15.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/17.0f * sinf(twopi * 17.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/19.0f * sinf(twopi * 19.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/21.0f * sinf(twopi * 21.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/23.0f * sinf(twopi * 23.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/25.0f * sinf(twopi * 25.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/27.0f * sinf(twopi * 27.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/29.0f * sinf(twopi * 29.0f *(float)i/size);
                *(x->wavetable+i) += 1.0f/31.0f * sinf(twopi * 31.0f *(float)i/size);
//            }
        }
        return (x);
    }
    
    // to-do
    // since we allocated some memory, we need a delete function
    static void pool_vox_free(t_pool_vox *x)
    {
        free(x->wavetable);
        
        //clear inlets
        inlet_free(x->x_in2); //used for sat
        inlet_free(x->x_in3); //used for value
        
        //clear outlets
        outlet_free(x->x_out1);     outlet_free(x->x_out2);
        outlet_free(x->x_out3);     outlet_free(x->x_out4);
        outlet_free(x->x_out5);     outlet_free(x->x_out6);
        outlet_free(x->x_out7);     outlet_free(x->x_out8);
        outlet_free(x->x_out9);     outlet_free(x->x_out10);
        outlet_free(x->x_out11);    outlet_free(x->x_out12);
        outlet_free(x->x_out13);    outlet_free(x->x_out14);
        outlet_free(x->x_out15);    outlet_free(x->x_out16);
    }
    
    /* this routine, which must have exactly this name (with the "~" replaced
     by "_tilde) is called when the code is first loaded, and tells Pd how
     to build the "class". */
    
    //class definition
    void pool_vox_tilde_setup(void)
    {
        
        // signal/float inlet 1 for hue, azimuth/elev = creatio args. 16 outputs.
        // creation args cannot be changed via inlets.
        pool_vox_class =        class_new(gensym("pool_vox~"),
                                          (t_newmethod)pool_vox_new, //creator class
                                          (t_method)pool_vox_free, //destructor class
                                          sizeof(t_pool_vox), //allocation
                                          CLASS_DEFAULT,//appearance
                                          A_DEFFLOAT,//creation arg (azimuth)
                                          A_DEFFLOAT,//creation arg (elevation)
                                          0);//termination
        
        /* here we tell Pd about the "dsp" method, which is called back
         when DSP is turned on. */
        class_addmethod(pool_vox_class, //specify class to add method to
                        (t_method)pool_vox_dsp, //cast to pd t_method
                        gensym("dsp"), //ID for dsp objects
                        A_CANT,//prevent crashes
                        0); //termination
        
        /* this is magic to declare that the leftmost, "main" inlet
         takes signals; other signal inlets are done differently... */
        CLASS_MAINSIGNALIN(pool_vox_class, t_pool_vox, x_f);
        //third argument of CLASS_MAINSIGNALIN is a (dummy-)floating point-variable of the data space, that is needed to automatically convert “float”-messages into signals if no signal is present at the signal-inlet.
        
    }
