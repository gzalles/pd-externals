#include "m_pd.h"
#include <stdlib.h>
#include <math.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

/* ------------------------ softclip~ ----------------------------- */

/* tilde object to do a simple soft clipping algorithm. */

static t_class *softclip_class;

typedef struct _softclip
{
    t_object x_obj; 	/* obligatory header */
    t_float gain;    	/* floats on the first inlet set the gain */
	t_float sample_rate;
    t_float freq;
} t_softclip;

    /* this is the actual performance routine which acts on the samples.
    It's called with a single pointer "w" which is our location in the
    DSP call list.  We return a new "w" which will point to the next item
    after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
    (no use to us), w[1] and w[2] are the input and output vector locations,
    and w[3] is the number of points to calculate. */
static t_int *softclip_perform(t_int *w)
{
	t_softclip *x = (t_softclip *)(w[1]);
    t_float *in = (t_float *)(w[2]);
    t_float *out = (t_float *)(w[3]);
    int n = (int)(w[4]);

	// i like counting from zero, so i use samplenumber to count the offset from
	// the start of the in and out blocks
	int samplenumber = 0;
	float a = -0.5f;
	float b = 0.0;
	float c = 1.5f;
	float d = 0.0;
	float ingain;

    while (n--)
    {
		ingain = x->gain * *(in+samplenumber);
		if(ingain > 1)
			*(out+samplenumber) = 1.0f;
		else if(ingain < -1)
			*(out+samplenumber) = -1.0f;
		else
			*(out+samplenumber) = a * ingain * ingain * ingain
							+ b * ingain * ingain
							+ c * ingain
							+ d;
		samplenumber++;
    }
    return (w+5);
}

    /* called to start DSP.  Here we call Pd back to add our perform
    routine to a linear callback list which Pd in turn calls to grind
    out the samples. */
static void softclip_dsp(t_softclip *x, t_signal **sp)
{
	x->sample_rate = sp[0]->s_sr;
    dsp_add(softclip_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

static void *softclip_new(void)
{
    t_softclip *x = (t_softclip *)pd_new(softclip_class);
    x->freq = 432.0f;
   outlet_new(&x->x_obj, gensym("signal"));
    return (x);
}

    /* this routine, which must have exactly this name (with the "~" replaced
    by "_tilde) is called when the code is first loaded, and tells Pd how
    to build the "class". */
void softclip_tilde_setup(void)
{
    softclip_class = class_new(gensym("softclip~"), (t_newmethod)softclip_new, 0,
    	sizeof(t_softclip), 0, A_DEFFLOAT, 0);
	
	/* this is magic to declare that the leftmost, "main" inlet
	    takes signals; other signal inlets are done differently... */
	/* also installs gain as the leftmost inlet float */
    CLASS_MAINSIGNALIN(softclip_class, t_softclip, gain);
	
	/* here we tell Pd about the "dsp" method, which is called back
	when DSP is turned on. */
    class_addmethod(softclip_class, (t_method)softclip_dsp, gensym("dsp"), (t_atomtype)0);
}