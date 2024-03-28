"""
-=-=-=-

Equations for calculating groups of diagrammatic contributions for single oscillator systems.

Vmax = 4
Tmax = 4
Minimum Coefficient: N/A

Generated 315 total diagrams in 0.023s.

-=-=-=-
"""


from scipy.special import factorial as fac

invf2 = 1/fac(2)
invf3 = 1/fac(3)
invf4 = 1/fac(4)

def get_diagrams(unpert, normal_ordered_coefficients, amplitudes, angu_freq, mode):

	contrib = dict({})

	hamil_delta, hamil_lambda, hamil_gamma = unpert
	R = 2*hamil_lambda
	Q_0_0, Q_0_1, Q_0_2, Q_0_3, Q_0_4 = normal_ordered_coefficients[0]
	Q_1_0, Q_1_1, Q_1_2, Q_1_3, _     = normal_ordered_coefficients[1]
	Q_2_0, Q_2_1, Q_2_2, _    , _     = normal_ordered_coefficients[2]
	Q_3_0, Q_3_1, _    , _    , _     = normal_ordered_coefficients[3]
	Q_4_0, _    , _    , _    , _     = normal_ordered_coefficients[4]

	t1, t2, t3, t4 = amplitudes

	if mode in ["all", "closed"]:

		#Ground State | 12 entries.
		contrib["0,0"] = (Q_0_0 + (1/2)*angu_freq + hamil_gamma) + \
			Q_0_1*1*t1 + \
			(Q_0_2 + R)*invf2*1*t2 + \
			(Q_0_2 + R)*invf2*1*t1*t1 + \
			Q_0_3*invf3*1*t3 + \
			Q_0_3*invf2*1*t1*t2 + \
			Q_0_3*invf3*1*t1*t1*t1 + \
			Q_0_4*invf4*1*t4 + \
			Q_0_4*invf3*1*t1*t3 + \
			Q_0_4*invf2*invf2*invf2*1*t2*t2 + \
			Q_0_4*invf2*invf2*1*t1*t1*t2 + \
			Q_0_4*invf4*1*t1*t1*t1*t1


	if mode in ["all", "open", "amplitude"]:

		#Amplitude t1 | 20 entries.
		contrib["1,0"] = Q_0_1*1*t2 + \
			Q_1_0 + \
			(Q_0_2 + R)*invf2*1*t3 + \
			(Q_0_2 + R)*1*t1*t2 + \
			(Q_1_1 + hamil_delta)*1*t1 + \
			Q_0_3*invf3*1*t4 + \
			Q_0_3*invf2*1*t1*t3 + \
			Q_0_3*invf2*invf2*1*t2*t2 + \
			Q_0_3*invf2*1*t1*t1*t2 + \
			Q_1_2*invf2*1*t2 + \
			Q_1_2*invf2*1*t1*t1 + \
			Q_0_4*invf3*1*t1*t4 + \
			Q_0_4*invf3*1*t2*t3 + \
			Q_0_4*invf2*invf2*1*t2*t3 + \
			Q_0_4*invf2*invf2*1*t1*t1*t3 + \
			Q_0_4*invf2*invf2*1*t1*t2*t2 + \
			Q_0_4*invf3*1*t1*t1*t1*t2 + \
			Q_1_3*invf3*1*t3 + \
			Q_1_3*invf2*1*t1*t2 + \
			Q_1_3*invf3*1*t1*t1*t1


		#Amplitude t2 | 30 entries.
		contrib["2,0"] = Q_0_1*0.5*t3 + \
			(Q_0_2 + R)*invf2*0.5*t4 + \
			(Q_0_2 + R)*0.5*t1*t3 + \
			(Q_0_2 + R)*invf2*1*t2*t2 + \
			(Q_1_1 + hamil_delta)*1*t2 + \
			invf2*(Q_2_0 + R) + \
			Q_0_3*invf2*0.5*t1*t4 + \
			Q_0_3*invf2*1*t2*t3 + \
			Q_0_3*invf2*0.5*t2*t3 + \
			Q_0_3*invf2*0.5*t1*t1*t3 + \
			Q_0_3*invf2*1*t1*t2*t2 + \
			Q_1_2*invf2*1*t3 + \
			Q_1_2*1*t1*t2 + \
			Q_2_1*0.5*t1 + \
			Q_0_4*invf3*1*t2*t4 + \
			Q_0_4*invf2*invf2*0.5*t2*t4 + \
			Q_0_4*invf2*invf2*invf2*1*t3*t3 + \
			Q_0_4*invf3*invf2*0.5*t3*t3 + \
			Q_0_4*invf2*invf2*0.5*t1*t1*t4 + \
			Q_0_4*invf2*1*t1*t2*t3 + \
			Q_0_4*invf2*0.5*t1*t2*t3 + \
			Q_0_4*invf2*invf3*1*t2*t2*t2 + \
			Q_0_4*invf3*0.5*t1*t1*t1*t3 + \
			Q_0_4*invf2*invf2*1*t1*t1*t2*t2 + \
			Q_1_3*invf3*1*t4 + \
			Q_1_3*invf2*1*t1*t3 + \
			Q_1_3*invf2*invf2*1*t2*t2 + \
			Q_1_3*invf2*1*t1*t1*t2 + \
			Q_2_2*invf2*0.5*t2 + \
			Q_2_2*invf2*0.5*t1*t1


		#Amplitude t3 | 34 entries.
		contrib["3,0"] = Q_0_1*0.16666666666666666*t4 + \
			(Q_0_2 + R)*0.16666666666666666*t1*t4 + \
			(Q_0_2 + R)*0.5*t2*t3 + \
			(Q_1_1 + hamil_delta)*0.5*t3 + \
			Q_0_3*invf2*0.5*t2*t4 + \
			Q_0_3*invf2*0.16666666666666666*t2*t4 + \
			Q_0_3*invf2*invf2*0.5*t3*t3 + \
			Q_0_3*invf2*0.16666666666666666*t1*t1*t4 + \
			Q_0_3*0.5*t1*t2*t3 + \
			Q_0_3*invf3*1*t2*t2*t2 + \
			Q_1_2*invf2*0.5*t4 + \
			Q_1_2*0.5*t1*t3 + \
			Q_1_2*invf2*1*t2*t2 + \
			Q_2_1*0.5*t2 + \
			invf3*Q_3_0 + \
			Q_0_4*invf3*0.5*t3*t4 + \
			Q_0_4*invf2*invf2*0.5*t3*t4 + \
			Q_0_4*invf3*0.16666666666666666*t3*t4 + \
			Q_0_4*invf2*0.5*t1*t2*t4 + \
			Q_0_4*invf2*0.16666666666666666*t1*t2*t4 + \
			Q_0_4*invf2*invf2*0.5*t1*t3*t3 + \
			Q_0_4*invf2*invf2*1*t2*t2*t3 + \
			Q_0_4*invf2*invf2*0.5*t2*t2*t3 + \
			Q_0_4*invf3*0.16666666666666666*t1*t1*t1*t4 + \
			Q_0_4*invf2*0.5*t1*t1*t2*t3 + \
			Q_0_4*invf3*1*t1*t2*t2*t2 + \
			Q_1_3*invf2*0.5*t1*t4 + \
			Q_1_3*invf2*1*t2*t3 + \
			Q_1_3*invf2*0.5*t2*t3 + \
			Q_1_3*invf2*0.5*t1*t1*t3 + \
			Q_1_3*invf2*1*t1*t2*t2 + \
			Q_2_2*invf2*0.5*t3 + \
			Q_2_2*0.5*t1*t2 + \
			Q_3_1*0.16666666666666666*t1


		#Amplitude t4 | 34 entries.
		contrib["4,0"] = (Q_0_2 + R)*0.16666666666666666*t2*t4 + \
			(Q_0_2 + R)*invf2*0.25*t3*t3 + \
			(Q_1_1 + hamil_delta)*0.16666666666666666*t4 + \
			Q_0_3*invf2*0.25*t3*t4 + \
			Q_0_3*invf2*0.16666666666666666*t3*t4 + \
			Q_0_3*0.16666666666666666*t1*t2*t4 + \
			Q_0_3*invf2*0.25*t1*t3*t3 + \
			Q_0_3*invf2*0.5*t2*t2*t3 + \
			Q_1_2*0.16666666666666666*t1*t4 + \
			Q_1_2*0.5*t2*t3 + \
			Q_2_1*0.25*t3 + \
			Q_0_4*invf2*invf2*invf2*0.25*t4*t4 + \
			Q_0_4*invf3*invf2*0.16666666666666666*t4*t4 + \
			Q_0_4*invf2*0.25*t1*t3*t4 + \
			Q_0_4*invf2*0.16666666666666666*t1*t3*t4 + \
			Q_0_4*invf2*invf2*0.5*t2*t2*t4 + \
			Q_0_4*invf2*invf2*0.16666666666666666*t2*t2*t4 + \
			Q_0_4*invf2*invf2*0.5*t2*t3*t3 + \
			Q_0_4*invf2*invf2*0.25*t2*t3*t3 + \
			Q_0_4*invf2*0.16666666666666666*t1*t1*t2*t4 + \
			Q_0_4*invf2*invf2*0.25*t1*t1*t3*t3 + \
			Q_0_4*invf2*0.5*t1*t2*t2*t3 + \
			Q_0_4*invf4*1*t2*t2*t2*t2 + \
			Q_1_3*invf2*0.5*t2*t4 + \
			Q_1_3*invf2*0.16666666666666666*t2*t4 + \
			Q_1_3*invf2*invf2*0.5*t3*t3 + \
			Q_1_3*invf2*0.16666666666666666*t1*t1*t4 + \
			Q_1_3*0.5*t1*t2*t3 + \
			Q_1_3*invf3*1*t2*t2*t2 + \
			Q_2_2*invf2*0.25*t4 + \
			Q_2_2*0.25*t1*t3 + \
			Q_2_2*invf2*0.5*t2*t2 + \
			Q_3_1*0.16666666666666666*t2 + \
			invf4*Q_4_0


	if mode in ["all", "open", "other"]:

		#=-=

		#Offset: -4 | 1 entries.
		contrib["0,4"] = invf4*Q_0_4

		#=-=

		#Offset: -3 | 2 entries.
		contrib["0,3"] = invf3*Q_0_3 + \
			invf3*Q_0_4*1*t1

		#=-=

		#Offset: -2 | 4 entries.
		contrib["0,2"] = invf2*(Q_0_2 + R) + \
			invf2*Q_0_3*1*t1 + \
			invf2*Q_0_4*invf2*1*t2 + \
			invf2*Q_0_4*invf2*1*t1*t1

		#Offset: -2 | 2 entries.
		contrib["1,3"] = invf3*Q_0_4*1*t2 + \
			invf3*Q_1_3

		#=-=

		#Offset: -1 | 7 entries.
		contrib["0,1"] = Q_0_1 + \
			(Q_0_2 + R)*1*t1 + \
			Q_0_3*invf2*1*t2 + \
			Q_0_3*invf2*1*t1*t1 + \
			Q_0_4*invf3*1*t3 + \
			Q_0_4*invf2*1*t1*t2 + \
			Q_0_4*invf3*1*t1*t1*t1

		#Offset: -1 | 5 entries.
		contrib["1,2"] = invf2*Q_0_3*1*t2 + \
			invf2*Q_1_2 + \
			invf2*Q_0_4*invf2*1*t3 + \
			invf2*Q_0_4*1*t1*t2 + \
			invf2*Q_1_3*1*t1

		#Offset: -1 | 1 entries.
		contrib["2,3"] = invf3*Q_0_4*0.5*t3

		#=-=

		#Offset: 0 | 11 entries.
		contrib["1,1"] = (Q_0_2 + R)*1*t2 + \
			(Q_1_1 + hamil_delta) + \
			Q_0_3*invf2*1*t3 + \
			Q_0_3*1*t1*t2 + \
			Q_1_2*1*t1 + \
			Q_0_4*invf3*1*t4 + \
			Q_0_4*invf2*1*t1*t3 + \
			Q_0_4*invf2*invf2*1*t2*t2 + \
			Q_0_4*invf2*1*t1*t1*t2 + \
			Q_1_3*invf2*1*t2 + \
			Q_1_3*invf2*1*t1*t1

		#Offset: 0 | 6 entries.
		contrib["2,2"] = invf2*Q_0_3*0.5*t3 + \
			invf2*Q_0_4*invf2*0.5*t4 + \
			invf2*Q_0_4*0.5*t1*t3 + \
			invf2*Q_0_4*invf2*1*t2*t2 + \
			invf2*Q_1_3*1*t2 + \
			invf2*invf2*Q_2_2

		#Offset: 0 | 1 entries.
		contrib["3,3"] = invf3*Q_0_4*0.16666666666666666*t4

		#=-=

		#Offset: 1 | 14 entries.
		contrib["2,1"] = (Q_0_2 + R)*0.5*t3 + \
			Q_0_3*invf2*0.5*t4 + \
			Q_0_3*0.5*t1*t3 + \
			Q_0_3*invf2*1*t2*t2 + \
			Q_1_2*1*t2 + \
			invf2*Q_2_1 + \
			Q_0_4*invf2*0.5*t1*t4 + \
			Q_0_4*invf2*1*t2*t3 + \
			Q_0_4*invf2*0.5*t2*t3 + \
			Q_0_4*invf2*0.5*t1*t1*t3 + \
			Q_0_4*invf2*1*t1*t2*t2 + \
			Q_1_3*invf2*1*t3 + \
			Q_1_3*1*t1*t2 + \
			Q_2_2*0.5*t1

		#Offset: 1 | 4 entries.
		contrib["3,2"] = invf2*Q_0_3*0.16666666666666666*t4 + \
			invf2*Q_0_4*0.16666666666666666*t1*t4 + \
			invf2*Q_0_4*0.5*t2*t3 + \
			invf2*Q_1_3*0.5*t3

		#=-=

		#Offset: 2 | 15 entries.
		contrib["3,1"] = (Q_0_2 + R)*0.16666666666666666*t4 + \
			Q_0_3*0.16666666666666666*t1*t4 + \
			Q_0_3*0.5*t2*t3 + \
			Q_1_2*0.5*t3 + \
			Q_0_4*invf2*0.5*t2*t4 + \
			Q_0_4*invf2*0.16666666666666666*t2*t4 + \
			Q_0_4*invf2*invf2*0.5*t3*t3 + \
			Q_0_4*invf2*0.16666666666666666*t1*t1*t4 + \
			Q_0_4*0.5*t1*t2*t3 + \
			Q_0_4*invf3*1*t2*t2*t2 + \
			Q_1_3*invf2*0.5*t4 + \
			Q_1_3*0.5*t1*t3 + \
			Q_1_3*invf2*1*t2*t2 + \
			Q_2_2*0.5*t2 + \
			invf3*Q_3_1

		#Offset: 2 | 3 entries.
		contrib["4,2"] = invf2*Q_0_4*0.16666666666666666*t2*t4 + \
			invf2*Q_0_4*invf2*0.25*t3*t3 + \
			invf2*Q_1_3*0.16666666666666666*t4

		#=-=

		#Offset: 3 | 11 entries.
		contrib["4,1"] = Q_0_3*0.16666666666666666*t2*t4 + \
			Q_0_3*invf2*0.25*t3*t3 + \
			Q_1_2*0.16666666666666666*t4 + \
			Q_0_4*invf2*0.25*t3*t4 + \
			Q_0_4*invf2*0.16666666666666666*t3*t4 + \
			Q_0_4*0.16666666666666666*t1*t2*t4 + \
			Q_0_4*invf2*0.25*t1*t3*t3 + \
			Q_0_4*invf2*0.5*t2*t2*t3 + \
			Q_1_3*0.16666666666666666*t1*t4 + \
			Q_1_3*0.5*t2*t3 + \
			Q_2_2*0.25*t3

		#Offset: 3 | 1 entries.
		contrib["5,2"] = invf2*Q_0_4*0.08333333333333333*t3*t4

		#=-=

		#Offset: 4 | 8 entries.
		contrib["5,1"] = Q_0_3*0.08333333333333333*t3*t4 + \
			Q_0_4*invf2*invf2*0.08333333333333333*t4*t4 + \
			Q_0_4*0.08333333333333333*t1*t3*t4 + \
			Q_0_4*invf2*0.16666666666666666*t2*t2*t4 + \
			Q_0_4*invf2*0.25*t2*t3*t3 + \
			Q_1_3*0.16666666666666666*t2*t4 + \
			Q_1_3*invf2*0.25*t3*t3 + \
			Q_2_2*0.08333333333333333*t4

		#Offset: 4 | 1 entries.
		contrib["6,2"] = invf2*Q_0_4*invf2*0.027777777777777776*t4*t4

		#=-=

		#Offset: 5 | 25 entries.
		contrib["5,0"] = (Q_0_2 + R)*0.08333333333333333*t3*t4 + \
			Q_0_3*invf2*invf2*0.08333333333333333*t4*t4 + \
			Q_0_3*0.08333333333333333*t1*t3*t4 + \
			Q_0_3*invf2*0.16666666666666666*t2*t2*t4 + \
			Q_0_3*invf2*0.25*t2*t3*t3 + \
			Q_1_2*0.16666666666666666*t2*t4 + \
			Q_1_2*invf2*0.25*t3*t3 + \
			Q_2_1*0.08333333333333333*t4 + \
			Q_0_4*invf2*invf2*0.08333333333333333*t1*t4*t4 + \
			Q_0_4*invf2*0.25*t2*t3*t4 + \
			Q_0_4*invf2*0.16666666666666666*t2*t3*t4 + \
			Q_0_4*invf2*0.08333333333333333*t2*t3*t4 + \
			Q_0_4*invf2*invf3*0.25*t3*t3*t3 + \
			Q_0_4*invf2*0.08333333333333333*t1*t1*t3*t4 + \
			Q_0_4*invf2*0.16666666666666666*t1*t2*t2*t4 + \
			Q_0_4*invf2*0.25*t1*t2*t3*t3 + \
			Q_0_4*invf3*0.5*t2*t2*t2*t3 + \
			Q_1_3*invf2*0.25*t3*t4 + \
			Q_1_3*invf2*0.16666666666666666*t3*t4 + \
			Q_1_3*0.16666666666666666*t1*t2*t4 + \
			Q_1_3*invf2*0.25*t1*t3*t3 + \
			Q_1_3*invf2*0.5*t2*t2*t3 + \
			Q_2_2*0.08333333333333333*t1*t4 + \
			Q_2_2*0.25*t2*t3 + \
			Q_3_1*0.08333333333333333*t3

		#Offset: 5 | 5 entries.
		contrib["6,1"] = Q_0_3*invf2*0.027777777777777776*t4*t4 + \
			Q_0_4*invf2*0.027777777777777776*t1*t4*t4 + \
			Q_0_4*0.08333333333333333*t2*t3*t4 + \
			Q_0_4*invf3*0.125*t3*t3*t3 + \
			Q_1_3*0.08333333333333333*t3*t4

		#=-=

		#Offset: 6 | 21 entries.
		contrib["6,0"] = (Q_0_2 + R)*invf2*0.027777777777777776*t4*t4 + \
			Q_0_3*invf2*0.027777777777777776*t1*t4*t4 + \
			Q_0_3*0.08333333333333333*t2*t3*t4 + \
			Q_0_3*invf3*0.125*t3*t3*t3 + \
			Q_1_2*0.08333333333333333*t3*t4 + \
			Q_0_4*invf2*invf2*0.08333333333333333*t2*t4*t4 + \
			Q_0_4*invf2*invf2*0.027777777777777776*t2*t4*t4 + \
			Q_0_4*invf2*invf2*0.125*t3*t3*t4 + \
			Q_0_4*invf2*invf2*0.08333333333333333*t3*t3*t4 + \
			Q_0_4*invf2*invf2*0.027777777777777776*t1*t1*t4*t4 + \
			Q_0_4*0.08333333333333333*t1*t2*t3*t4 + \
			Q_0_4*invf3*0.125*t1*t3*t3*t3 + \
			Q_0_4*invf3*0.16666666666666666*t2*t2*t2*t4 + \
			Q_0_4*invf2*invf2*0.25*t2*t2*t3*t3 + \
			Q_1_3*invf2*invf2*0.08333333333333333*t4*t4 + \
			Q_1_3*0.08333333333333333*t1*t3*t4 + \
			Q_1_3*invf2*0.16666666666666666*t2*t2*t4 + \
			Q_1_3*invf2*0.25*t2*t3*t3 + \
			Q_2_2*0.08333333333333333*t2*t4 + \
			Q_2_2*invf2*0.125*t3*t3 + \
			Q_3_1*0.027777777777777776*t4

		#Offset: 6 | 3 entries.
		contrib["7,1"] = Q_0_4*invf2*0.027777777777777776*t2*t4*t4 + \
			Q_0_4*invf2*0.041666666666666664*t3*t3*t4 + \
			Q_1_3*invf2*0.027777777777777776*t4*t4

		#=-=

		#Offset: 7 | 13 entries.
		contrib["7,0"] = Q_0_3*invf2*0.027777777777777776*t2*t4*t4 + \
			Q_0_3*invf2*0.041666666666666664*t3*t3*t4 + \
			Q_1_2*invf2*0.027777777777777776*t4*t4 + \
			Q_0_4*invf2*invf2*0.041666666666666664*t3*t4*t4 + \
			Q_0_4*invf2*invf2*0.027777777777777776*t3*t4*t4 + \
			Q_0_4*invf2*0.027777777777777776*t1*t2*t4*t4 + \
			Q_0_4*invf2*0.041666666666666664*t1*t3*t3*t4 + \
			Q_0_4*invf2*0.08333333333333333*t2*t2*t3*t4 + \
			Q_0_4*invf3*0.125*t2*t3*t3*t3 + \
			Q_1_3*invf2*0.027777777777777776*t1*t4*t4 + \
			Q_1_3*0.08333333333333333*t2*t3*t4 + \
			Q_1_3*invf3*0.125*t3*t3*t3 + \
			Q_2_2*0.041666666666666664*t3*t4

		#Offset: 7 | 1 entries.
		contrib["8,1"] = Q_0_4*invf2*0.013888888888888888*t3*t4*t4

		#=-=

		#Offset: 8 | 9 entries.
		contrib["8,0"] = Q_0_3*invf2*0.013888888888888888*t3*t4*t4 + \
			Q_0_4*invf2*invf3*0.013888888888888888*t4*t4*t4 + \
			Q_0_4*invf2*0.013888888888888888*t1*t3*t4*t4 + \
			Q_0_4*invf2*invf2*0.027777777777777776*t2*t2*t4*t4 + \
			Q_0_4*invf2*0.041666666666666664*t2*t3*t3*t4 + \
			Q_0_4*invf4*0.0625*t3*t3*t3*t3 + \
			Q_1_3*invf2*0.027777777777777776*t2*t4*t4 + \
			Q_1_3*invf2*0.041666666666666664*t3*t3*t4 + \
			Q_2_2*invf2*0.013888888888888888*t4*t4

		#Offset: 8 | 1 entries.
		contrib["9,1"] = Q_0_4*invf3*0.004629629629629629*t4*t4*t4

		#=-=

		#Offset: 9 | 5 entries.
		contrib["9,0"] = Q_0_3*invf3*0.004629629629629629*t4*t4*t4 + \
			Q_0_4*invf3*0.004629629629629629*t1*t4*t4*t4 + \
			Q_0_4*invf2*0.013888888888888888*t2*t3*t4*t4 + \
			Q_0_4*invf3*0.020833333333333332*t3*t3*t3*t4 + \
			Q_1_3*invf2*0.013888888888888888*t3*t4*t4

		#=-=

		#Offset: 10 | 3 entries.
		contrib["10,0"] = Q_0_4*invf3*0.004629629629629629*t2*t4*t4*t4 + \
			Q_0_4*invf2*invf2*0.006944444444444444*t3*t3*t4*t4 + \
			Q_1_3*invf3*0.004629629629629629*t4*t4*t4

		#=-=

		#Offset: 11 | 1 entries.
		contrib["11,0"] = Q_0_4*invf3*0.0023148148148148147*t3*t4*t4*t4

		#=-=

		#Offset: 12 | 1 entries.
		contrib["12,0"] = Q_0_4*invf4*0.0007716049382716049*t4*t4*t4*t4

	return contrib
