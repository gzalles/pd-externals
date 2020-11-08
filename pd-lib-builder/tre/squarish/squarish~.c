#include "m_pd.h"
#include <stdlib.h>
#include <math.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

/* ------------------------ squarish~ ----------------------------- */

#define WAVETABLESIZE 1024

static t_class *squarish_class;

typedef struct _squarish
{
    t_object x_obj; 	/* obligatory header */
    t_float x_f;    	/* place to hold inlet's value if it's set by message */
	t_float *wavetable;	/* a place to hold the squared values */
	t_float phase;
	t_float samplerate;
} t_squarish;

// some nice little interpolation routines
static inline float no_interpolate(t_squarish *x)
{
	int x_1 = x->phase * WAVETABLESIZE;

	return(x->wavetable[x_1 % WAVETABLESIZE]);
}
static inline float lin_interpolate(t_squarish *x)
{
	int x_1 = x->phase * WAVETABLESIZE;
	float y_1 = x->wavetable[x_1 % WAVETABLESIZE];
	float y_2 = x->wavetable[(x_1 + 1) % WAVETABLESIZE];

	return (y_2 - y_1) * ((x->phase * WAVETABLESIZE) - x_1) + y_1;
}

static inline float quad_interpolate(t_squarish *x)
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

static t_int *squarish_perform(t_int *w)
{
	t_squarish *x = (t_squarish *)(w[1]);
    t_float *freq = (t_float *)(w[2]);
    t_float *out = (t_float *)(w[3]);
    int n = (int)(w[4]);

	// i like counting from zero
	int blocksize = n;
	int i, sample = 0;
	float phaseincrement;
	float findex;
    float wphase;
	int	iindex;
	
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
    	*(out+sample) = x->phase;
		// now grab the sample from the table
		findex = WAVETABLESIZE * wphase;
		iindex = (int)findex;
//		*(out+sample) = *(x->wavetable + iindex);
//		*(out+sample) = lin_interpolate(x);
		*(out+sample) = quad_interpolate(x);
		sample++;
    }
    return (w+5);
}

    /* called to start DSP.  Here we call Pd back to add our perform
    routine to a linear callback list which Pd in turn calls to grind
    out the samples. */

static void squarish_dsp(t_squarish *x, t_signal **sp)
{
	// we'll initialize samplerate when starting up
	x->samplerate = sp[0]->s_sr;
    dsp_add(squarish_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

//constructor method
static void *squarish_new(void)
{
	float twopi, size;
	int i;
    
    //allocate memory for the object itself
    t_squarish *x = (t_squarish *)pd_new(squarish_class);
    
    outlet_new(&x->x_obj, gensym("signal"));
    
	// initialize variables
    x->x_f = 0.0f;
	x->phase = 0.0f;
	twopi = 8.0f * atanf(1.0f);
	size = (float)WAVETABLESIZE;
	
	// space for WAVETABLESIZE samples
	x->wavetable = (t_float *)malloc(WAVETABLESIZE * sizeof(t_float));
	
	// fill it up with a sine wave
	for(i = 0; i < WAVETABLESIZE; i++)
    {
        *(x->wavetable+i) = sinf(twopi * (float)i/size);
        *(x->wavetable+i) += 1.0f/3.0f * sinf(twopi * 3.0f *(float)i/size);
        *(x->wavetable+i) += 1.0f/5.0f * sinf(twopi * 5.0f *(float)i/size);
        *(x->wavetable+i) += 1.0f/7.0f * sinf(twopi * 7.0f *(float)i/size);
        *(x->wavetable+i) += 1.0f/9.0f * sinf(twopi * 9.0f *(float)i/size);
    }
    
    return (x);
}

// since we allocated some memory, we need a delete function
static void squarish_free(t_squarish *x)
{
	free(x->wavetable);
}

    /* this routine, which must have exactly this name (with the "~" replaced
    by "_tilde) is called when the code is first loaded, and tells Pd how
    to build the "class". */

//class definition
void squarish_tilde_setup(void)
{
    squarish_class = class_new(gensym("squarish~"), (t_newmethod)squarish_new, (t_method)squarish_free,
    	sizeof(t_squarish), 0, A_DEFFLOAT, 0);
	    /* this is magic to declare that the leftmost, "main" inlet
	    takes signals; other signal inlets are done differently... */
    CLASS_MAINSIGNALIN(squarish_class, t_squarish, x_f);
    	/* here we tell Pd about the "dsp" method, which is called back
	when DSP is turned on. */
    class_addmethod(squarish_class, (t_method)squarish_dsp, gensym("dsp"), 0);
}
