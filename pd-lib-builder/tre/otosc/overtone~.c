#include "m_pd.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

void *overtone_tilde_class;

#define WAVETABLESIZE 8192
#define HARMONICS 120

typedef struct _overtone_tilde
{
    t_object x_obj; 	/* obligatory header */
	t_float bands;
	t_float gain;
	t_float slope;
	t_float *wavetable;
	t_float *position;
	t_float *amplitude;
	t_float *slopeamp;
	t_float *freqoffset;
	t_float *phase;
	t_float frequency;
	t_float srate;
	t_float nyquist;
	t_float oneoversrate;
	t_outlet *outone;
	t_outlet *outtwo;

} t_overtone_tilde;

void *overtone_tilde_new(t_symbol *s, short argc, t_atom *argv);
void overtone_tilde_free(t_overtone_tilde *x);
t_int *overtone_tilde_perform(t_int *w);
void overtone_tilde_dsp(t_overtone_tilde *x, t_signal **sp);
void overtone_tilde_setband(t_overtone_tilde *x, float f);
void overtone_tilde_setslope(t_overtone_tilde *x, float f);
void overtone_tilde_setfrequency(t_overtone_tilde *x, float f);
void overtone_tilde_assist(t_overtone_tilde *x, void *b, long m, long a, char *s);
void overtone_tilde_bandset(t_overtone_tilde *x, t_symbol *message, short argc, t_atom *argv); 
inline float no_interpolate(t_overtone_tilde *x, int band);


void overtone_tilde_setup(void)
{
    overtone_tilde_class = class_new(gensym("overtone~"), (t_newmethod)overtone_tilde_new, (t_method)overtone_tilde_free,
    	sizeof(t_overtone_tilde), 0, A_GIMME, 0);
    class_addfloat(overtone_tilde_class, (t_method)overtone_tilde_setfrequency);
    	/* here we tell Pd about the "dsp" method, which is called back
	when DSP is turned on. */
    class_addmethod(overtone_tilde_class, (t_method)overtone_tilde_dsp, gensym("dsp"), 0);
    class_addmethod(overtone_tilde_class, (t_method)overtone_tilde_bandset, gensym("set"), A_GIMME, 0);
    class_addmethod(overtone_tilde_class, (t_method)overtone_tilde_setband, gensym("band"), A_FLOAT, 0);
    class_addmethod(overtone_tilde_class, (t_method)overtone_tilde_setslope, gensym("slope"), A_FLOAT, 0);
}


void *overtone_tilde_new(t_symbol *s, short argc, t_atom *argv)
{
	long i;
	
	float oneoversize = 1.0f/WAVETABLESIZE;
	float twoPi = 8.0f * atan(1.0f);

    t_overtone_tilde *x = (t_overtone_tilde *)pd_new(overtone_tilde_class);

	x->outone = outlet_new(&x->x_obj, gensym("signal"));
	x->outtwo = outlet_new(&x->x_obj, gensym("signal"));
    
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("band"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("slope"));
    
	x->wavetable = (float *)malloc(sizeof(float) * WAVETABLESIZE);
	for(i=0; i<WAVETABLESIZE; i++)
		*(x->wavetable + i) = sin(twoPi * (float)i*oneoversize);
		
	x->position = (float *)malloc(sizeof(float) * HARMONICS);
	for(i=0; i<HARMONICS; i++)
		*(x->position +i) = ((float)random())/2147483648.0f + 0.5f;
		
	x->phase = (float *)malloc(sizeof(float) * HARMONICS);
	for(i=0; i<HARMONICS; i++)
//		x->phase[i] = ((float)random()*8192.0f)/2147483648.0f;
		x->phase[i] = 0.0f;
		
	x->amplitude = (float *)malloc(sizeof(float) * HARMONICS);
	x->slopeamp = (float *)malloc(sizeof(float) * HARMONICS);
	x->freqoffset = (float *)malloc(sizeof(float) * HARMONICS);
	
    for(i=0; i<HARMONICS; i++)
		x->amplitude[i] = x->slopeamp[i] = x->freqoffset[i] = 1.0f;
		
	x->bands = 2.0f;
    if(x->bands > 1.0f)
    	x->gain = 1.0f/(x->bands);
    else
    	x->gain = 1.0f;

	post("overtone_tilde: oscillator with up to 120 partials");
    return (x);
}

// free memory here
void overtone_tilde_free(t_overtone_tilde *x)
{
	free(x->wavetable);
	free(x->position);
	free(x->phase);
	free(x->amplitude);
	free(x->slopeamp);
	free(x->freqoffset);
}

void overtone_tilde_dsp(t_overtone_tilde *x, t_signal **sp)
{
	x->srate =  sp[0]->s_sr;
	x->nyquist = x->srate * 0.5f;
	x->oneoversrate = 1.0/x->srate;
	dsp_add(overtone_tilde_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

void overtone_tilde_bandset(t_overtone_tilde *x, t_symbol *message, short argc, t_atom *argv)
{
	long band;
	t_float amount;
	
	if(argc != 3)
	{
		post("set requires 3 arguments: freq|pan|amp (name), bandnumber (float), amount (float)");
		return;
	}
	// first argument
	if(argv[0].a_type == A_SYMBOL)
	{
		if(strcmp(argv[0].a_w.w_symbol->s_name, "freq") == 0)
		{
			if(argv[1].a_type == A_FLOAT)
			{
				if(argv[2].a_type == A_FLOAT)
				{
					band = argv[1].a_w.w_float;
					if(band < 0)
						band = 0;
					if(band > 119)
						band = 119;
					amount = argv[2].a_w.w_float;
					if(amount < 0.01f)
						amount = 0.01f;
					if(amount > x->nyquist)
						amount = x->nyquist;
					x->freqoffset[band] = amount / (x->frequency * (band + 1));
				}
				else
				{
					post("the third argument to \"set freq\" should be a float");
					return;
				}
			}
			else
			{
				post("the second argument to \"set freq\" should be an float");
				return;
			}
		}
		else if(strcmp(argv[0].a_w.w_symbol->s_name, "pan") == 0)
		{
			if(argv[1].a_type == A_FLOAT)
			{
				if(argv[2].a_type == A_FLOAT)
				{
					band = argv[1].a_w.w_float;
					if(band < 0)
						band = 0;
					if(band > 119)
						band = 119;
					amount = argv[2].a_w.w_float;
					if(amount < 0.0f)
						amount = 0.0f;
					if(amount > 1.0f)
						amount = 1.0f;
					x->position[band] = amount;
				}
				else
				{
					post("the third argument to \"set pan\" should be a float");
					return;
				}
			}
			else
			{
				post("the second argument to \"set pan\" should be an float");
				return;
			}
		}
		else if(strcmp(argv[0].a_w.w_symbol->s_name, "amp") == 0)
		{
			if(argv[1].a_type == A_FLOAT)
			{
				if(argv[2].a_type == A_FLOAT)
				{
					band = argv[1].a_w.w_float;
					if(band < 0)
						band = 0;
					if(band > 119)
						band = 119;
					amount = argv[2].a_w.w_float;
					if(amount < 0.0f)
						amount = 0.0f;
					if(amount > 1.0f)
						amount = 1.0f;
					x->amplitude[band] = amount;
				}
				else
				{
					post("the third argument to \"set amp\" should be a float");
					return;
				}
			}
			else
			{
				post("the second argument to \"set amp\" should be an float");
				return;
			}
		}
		else
		{
			post("the first argument to \"set\" should be \"freq\", \"pan\", or \"amp\".");
			return;
		}
	}
	else
	{
		post("the first argument to \"set\" should be \"freq\", \"pan\", or \"amp\".");
		return;
	}

}

void overtone_tilde_setband(t_overtone_tilde *x, float f)
{
    x->bands = f;
    if(x->bands > 119.0f)
    	x->bands = 119.0f;
    if(x->bands < 0.0f)
    	x->bands = 0.0f;
    if(x->bands > 1.0f)
    	x->gain = 1.0f/(x->bands);
    else
    	x->gain = 1.0f;
}

void overtone_tilde_setfrequency(t_overtone_tilde *x, float f)
{
    x->frequency = f;
    if(x->frequency > x->nyquist)
    	x->frequency = x->nyquist;
    if(x->frequency < 0.0f)
    	x->frequency = 0.0f;
}

void overtone_tilde_setslope(t_overtone_tilde *x, float f)
{
	float slopeBasis, sum, adjust;
	long band;
	
    x->slope = (f * 2.0f) - 1.0f;
    if(x->slope > 1.0f)
    	x->slope = 1.0f;
    if(x->slope < -1.0f)
    	x->slope = -1.0f;

	slopeBasis = log10((float)HARMONICS) * 20.0f * 0.5f;
	for(band = 0; band < HARMONICS; band++)
		x->slopeamp[band] = pow(10.0f, ((log10((float)(band + 1.0f)) * 20.0f) - slopeBasis) * (x->slope/slopeBasis));
		sum = 0.0f;
	for(band = 0; band < HARMONICS; band++)
		sum += x->slopeamp[band];
	adjust = (float)HARMONICS/sum;
	for(band = 0; band < HARMONICS; band++)
		x->slopeamp[band] = x->slopeamp[band] * adjust;
	
}

inline float no_interpolate(t_overtone_tilde *x, int band)
{
	int x_1 = x->phase[band];

	return(x->wavetable[x_1 % WAVETABLESIZE]);
}

t_int *overtone_tilde_perform(t_int *w)
{
	t_overtone_tilde *x = (t_overtone_tilde *) (w[1]);
	float *outupper = (t_float *)(w[2]);
	float *outlower = (t_float *)(w[3]);
	int n = (int)(w[4]);
	
	long i, band;
	float phaseinc;
	float sample, scaleupper, scalelower, scale;
	long longband = (long)(x->bands);
	float fracband = x->bands - (float)longband;
	float gain = x->gain;
	float tablefreq = x->frequency * WAVETABLESIZE *  x->oneoversrate;

	for (i = 0; i < n; i++)
	{
		*(outlower+i) = 0;
		*(outupper+i) = 0;
	}
	// low bands
	for(band = 0; band < longband ; band++)
	{
		while (x->phase[band] >= WAVETABLESIZE)
			x->phase[band] -= WAVETABLESIZE;
		while (x->phase[band] < 0)
			x->phase[band] += WAVETABLESIZE;
		phaseinc = (band + 1.0f) * tablefreq * x->freqoffset[band];
		scale = gain * x->amplitude[band] * x->slopeamp[band];
		scalelower = scale * x->position[band];
		scaleupper = scale * (1.0f - x->position[band]);
		for (i = 0; i < n; i++)
		{
			x->phase[band] += phaseinc;
			sample = no_interpolate(x, band);
			*(outlower+i) += sample * scalelower;
			*(outupper+i) += sample * scaleupper;
		}
	}
	// high band
	while (x->phase[longband] >= WAVETABLESIZE)
		x->phase[longband] -= WAVETABLESIZE;
	while (x->phase[longband] < 0)
		x->phase[longband] += WAVETABLESIZE;
	phaseinc = (longband + 1.0f) * tablefreq * x->freqoffset[band];
	scale = gain * x->amplitude[longband] * x->slopeamp[longband] * fracband;
	scalelower = scale * x->position[longband];
	scaleupper = scale * (1.0f - x->position[longband]);
	for (i = 0; i < n; i++)
	{
		x->phase[longband] += phaseinc;
		sample = no_interpolate(x, longband);
		*(outlower+i) += sample * scalelower;
		*(outupper+i) += sample * scaleupper;
	}
	return (w+5);
}

