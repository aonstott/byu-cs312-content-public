import random

# This is main function that is connected to the Test button. You don't need to touch it.
def prime_test(N, k):
	return fermat(N,k), miller_rabin(N,k)

# You will need to implement this function and change the return value.
def mod_exp(x, y, N): 
	if y == 0:
		return 1
	z = mod_exp(x, y//2, N)

	if y % 2 == 0:
		return z**2 % N
	else:
		return x * z**2 % N


# You will need to implement this function and change the return value.   
def fprobability(k):

    return 1 - (1/2)**k

# You will need to implement this function and change the return value.   
def mprobability(k):
	return 1 - (1/4)**k
    

# You will need to implement this function and change the return value, which should be
# either 'prime' or 'composite'.
#
# To generate random values for a, you will most likley want to use
# random.randint(low,hi) which gives a random integer between low and
# hi, inclusive.
def fermat(N,k):
	a_vals = []
	for i in range (0, k):
		a_vals.append(random.randint(2, N-1))
	for a in a_vals:
		if mod_exp(a, N-1, N) != 1:
			return 'composite'
		
	return 'prime'

	

# You will need to implement this function and change the return value, which should be
# either 'prime' or 'composite'.
#
# To generate random values for a, you will most likley want to use
# random.randint(low,hi) which gives a random integer between low and
#  hi, inclusive.
def miller_rabin(N,k):
	a_vals = [random.randint(2, N - 1) for _ in range(k)]

	for a in a_vals:
		exp = N - 1

		if N % 2 == 0:
			return 'composite'
		while exp % 2 == 0:
			result = mod_exp(a, exp, N)
			if result == N - 1:
				break
			if result != 1:
				return 'composite'
			exp //= 2

	return 'prime'

#test = fermat(5, 3)
#main function to test in terminal 
def main():
	print(fermat(5, 3))
	print(miller_rabin(156, 3))
	print(mod_exp(2, 3, 5))
	print(fprobability(5))
	print(mprobability(5))

if __name__ == "__main__":
	main()
