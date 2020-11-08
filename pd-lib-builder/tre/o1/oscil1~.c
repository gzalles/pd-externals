#include "m_pd.h"
#include <stdlib.h>
#include <math.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

/* ------------------------ oscil1~ ----------------------------- */

#define WAVETABLESIZE 1024 //here is the macro def for wavetable size

static t_class *oscil1_class;

typedef struct _oscil1
{
    t_object x_obj; 	/* obligatory header */
    t_float x_f;    	/* place to hold inlet's value if it's set by message */
	t_float *wavetable;	/* a place to hold the squared values */
	t_float phase;
	t_float samplerate;
    
} t_oscil1;

// some nice little interpolation routines
static inline float no_interpolate(t_oscil1 *x)
{
	int x_1 = x->phase * WAVETABLESIZE;

	return(x->wavetable[x_1 % WAVETABLESIZE]);
}

static inline float lin_interpolate(t_oscil1 *x, float harm)
{
	int x_1 = x->phase * WAVETABLESIZE * harm;
	float y_1 = x->wavetable[x_1 % WAVETABLESIZE];
	float y_2 = x->wavetable[(x_1 + 1) % WAVETABLESIZE];

	return (y_2 - y_1) * ((x->phase * WAVETABLESIZE * harm) - x_1) + y_1;
}

static inline float quad_interpolate(t_oscil1 *x)
{
	int truncphase = (int) (x->phase * WAVETABLESIZE);
	float fr = (x->phase * WAVETABLESIZE) - ((float) truncphase);
	float inm1 = x->wavetable[(truncphase - 1) % WAVETABLESIZE];
	float in   = x->wavetable[(truncphase + 0) % WAVETABLESIZE];
	float inp1 = x->wavetable[(truncphase + 1) % WAVETABLESIZE];
	float inp2 = x->wavetable[(truncphase + 2) % WAVETABLESIZE];

	return in + 0.5 * fr * (inp1 - inm1 +
	 fr * (4.0 * inp1 + 2.0 * inm1 - 5.0 * in - inp2 +
	 fr * (3.0 * (in - inp1) - inm1 + inp2)));
}

    /* this is the actual performance routine which acts on the samples.
    It's called with a single pointer "w" which is our location in the
    DSP call list.  We return a new "w" which will point to the next item
    after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
    (no use to us), w[1] and w[2] are the input and output vector locations,
    and w[3] is the number of points to calculate. */
static t_int *oscil1_perform(t_int *w)
{
	t_oscil1 *x = (t_oscil1 *)(w[1]);   //pointer to the object itself
    t_float *freq = (t_float *)(w[2]);  //value being passed by sig~ object
    t_float *out = (t_float *)(w[3]);   //the output
    int n = (int)(w[4]);                //block size

	// i like counting from zero
	int blocksize = n;
	int i, j, sample = 0;//init
    
	float phaseincrement, triflip;
	
    float findex; //float index
	int	iindex; //integer index
	
    while (n--)
    {
		// first we need to calculate the phase increment from the frequency
		// and sample rate - this is the number of cycles per sample
		// freq = cyc/sec, sr = samp/sec, phaseinc = cyc/samp = freq/sr
		phaseincrement = *(freq+sample)/x->samplerate;
		
		// now, increment the phase and make sure it doesn't go over 1.0
		x->phase = x->phase + phaseincrement;
		
        //if the phase is greater than or equal to one subtract one
        while(x->phase >= 1.0f)
			x->phase = x->phase - 1.0f;
        //if is less than zero add 1, this is in case user puts negative values in sig~ obj.
		while(x->phase < 0.0f)
			x->phase = x->phase + 1.0f;
        
		// now grab the sample from the table using the float index
		findex = WAVETABLESIZE * x->phase;
		
        //round to an integer, casting, this is why we interpolate...
        iindex = (int)findex;
        
        *(out+sample) = 0.0f;
        // a not so efficient triangle wave generator. first 9 harmonics
        // invert every other harmonic to create the triangle shape
        triflip = 1.0f;
        
        for(j = 1; j < 17; j += 2)
        {
            *(out+sample) += ((1/j^2) * triflip * lin_interpolate(x, (float)j)/(float)(j*j));
            triflip *= -1.0f;
        }
        
		sample++;
    }
    return (w+5);
}

    /* called to start DSP.  Here we call Pd back to add our perform
    routine to a linear callback list which Pd in turn calls to grind
    out the samples. */
static void oscil1_dsp(t_oscil1 *x, t_signal **sp)
{
	// we'll initialize samplerate when starting up sp[0] = sample rate
	x->samplerate = sp[0]->s_sr;
    //pass 4 params
    dsp_add(oscil1_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

//creation
static void *oscil1_new(void)
{
	float twopi, size;
	int i;
    
    t_oscil1 *x = (t_oscil1 *)pd_new(oscil1_class);
    
    outlet_new(&x->x_obj, gensym("signal"));
	
    // initialize variables
    x->x_f = 0.0f;
	x->phase = 0.0f;
	twopi = 8.0f * atanf(1.0f);//trick to get 2 pi easily
	size = (float)WAVETABLESIZE; //macro defined at beginning, casted to float
	
	// space for WAVETABLESIZE samples, memory allocation
	x->wavetable = (t_float *)malloc(WAVETABLESIZE * sizeof(t_float));
	
	// fill it up with a sine wave
	for(i = 0; i < WAVETABLESIZE; i++)
        *(x->wavetable+i) = sinf(twopi * (float)i/size);//using math.h
    return (x);
}

// since we allocated some memory, we need a delete function
static void oscil1_free(t_oscil1 *x)
{
	free(x->wavetable);
}

/* this routine, which must have exactly this name (with the "~" replaced
by "_tilde) is called when the code is first loaded, and tells Pd how to build the "class". */
void oscil1_tilde_setup(void)
{
    oscil1_class = class_new(gensym("oscil1~"), (t_newmethod)oscil1_new, (t_method)oscil1_free,
    	sizeof(t_oscil1), 0, A_DEFFLOAT, 0);
    
    /* this is magic to declare that the leftmost, "main" inlet
    takes signals; other signal inlets are done differently... */
    CLASS_MAINSIGNALIN(oscil1_class, t_oscil1, x_f);
    
    /* here we tell Pd about the "dsp" method, which is called back
	when DSP is turned on. */
    class_addmethod(oscil1_class, (t_method)oscil1_dsp, gensym("dsp"), 0);
}
