#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cfloat>

# define M_PI           3.14159265358979323846

using std::cout;
using std::endl;
using std::setw;
using std::vector;
float xmin = 0;
float xmax = M_PI;

size_t K = 100;
std::vector<float> getFtilda() {
	vector<float> Ftilda(K);
	for (size_t k = 0; k < K; ++k) {
		Ftilda[k] = sin(xmin + k * (xmax - xmin) / K) + 0.5 + (float)(rand()) / RAND_MAX * 0.5 - 0.25;
	}
	return Ftilda;
}
std::vector<float> getAlpha3() {
	std::vector<float> alpha(3);
	alpha[1] = (float)(rand()) / RAND_MAX;
	alpha[0] = alpha[2] = (1 - alpha[1]) / 2;
	return alpha;
}
std::vector<float> getAlpha5() {
	std::vector<float> alpha(5);
	alpha[2] = (float)(rand()) / RAND_MAX;
	float random = (float)(rand()) / RAND_MAX;
	while (random > 1 - alpha[2]) {
		random = (float)(rand()) / RAND_MAX;
	}
	alpha[1] = alpha[3] = random / 2;
	alpha[0] = alpha[4] = (1 - alpha[2] - 2 * alpha[1]) / 2;
	return alpha;
}
std::vector<float> getFdash(vector<float> Ftilda, vector<float> alpha) {
	std::vector<float> Fdash(K);
	size_t M = (alpha.size() - 1) / 2;
	for (size_t k = M; k < K - M; k++) {
		Fdash[k] = 0;
		for (size_t j = 0; j < alpha.size(); j++) {
			Fdash[k] += Ftilda[j] * Ftilda[k + j - M] * alpha[j];
		}
	}
	return Fdash;
}
float getW(vector<float> Fdash) {
	float w = 0;
	for (size_t i = 1; i < K; ++i) {
		w += abs(Fdash[i] - Fdash[i - 1]);
	}
	return w;
}
float getDelta(vector<float> Fdash, vector<float> Ftilda) {
	float delta = 0;
	for (size_t i = 0; i < K; ++i) {
		delta += abs(Fdash[i] - Ftilda[i]);
	}
	delta /= K;
	return delta;
}
void print(float lyambda, float J, vector<float> alpha, float w, float delta) {
	if (lyambda == 0) {
		if (alpha.size() == 3) {
			cout
				<< "+---------+---------------+---------------+---------------------------------------------+---------------+---------------+"
				<< endl;
			cout
				<< "| lyambda | Distance | J | alpha | w | delta |"
				<< endl;
			cout
				<< "+---------+---------------+---------------+---------------------------------------------+---------------+---------------+"
				<< endl;
		}
		else {
			cout
				<< "+---------+---------------+---------------+---------------------------------------------------------------------------+--------------- +-------------- - +"<< endl;
			cout
				<< "| lyambda | Distance | J | alpha | w | delta | " << endl;
			cout << "+---------+---------------+---------------+---------------------------------------------------------------------------+--------------- +-------------- - +"<< endl;
		}
	}

	cout << "|" << setw(9) << lyambda << "|" << setw(15) << abs(w) + abs(delta) << "|" << setw(15) << J << "|";
	for (int i = 0; i < alpha.size(); ++i) {
		cout << setw(15) << alpha[i];
	}
	cout << "|" << setw(15) << w << "|" << setw(15) << delta << "|" << endl;
	if (lyambda > 1) {
		if (alpha.size() == 3) {
			cout
				<< "+---------+---------------+---------------+---------------------------------------------+---------------+---------------+"
				<< endl;
		}
		else {
			cout
				<< "+---------+---------------+---------------+---------------------------------------------------------------------------+--------------- +-------------- - +"
				<< endl;
		}
		cout << endl;
	}
}
void print_best(float lyambda, float J, float w, float delta) {
	cout << "+---------+---------------+---------------+---------------+" << endl;
	cout << "| lyambda | J | w | delta |" << endl;
	cout << "+---------+---------------+---------------+---------------+" << endl;
	cout << "|" << setw(9) << lyambda << "|" << setw(15) << J << "|" << setw(15) << w << "|" << setw(15) << delta <<
		"|"
		<< endl;
	cout << "+---------+---------------+---------------+---------------+" << endl << endl;
}
int main() {
	srand(time(NULL));
	vector<float> Ftilda(getFtilda());
	int r = 3;
	for (size_t z = 0; z < 2; ++z) {
		float best_lyambda = FLT_MAX;
		float best_lyambda_J = FLT_MAX;
		float best_lyambda_w = FLT_MAX;
		float best_lyambda_delta = FLT_MAX;
		vector<float> best_alpha;
		for (float lyambda = 0; lyambda < 1.1; lyambda += 0.1) {
			bool flag = 0;
			float J = FLT_MAX;
			float w = FLT_MAX;
			float delta = FLT_MAX;
			vector<float> alpha;
			for (size_t i = 0; i < 10; ++i) {
				vector<float> buff_alpha;
				if (r == 3) {
					buff_alpha = getAlpha3();
				}
				else {
					buff_alpha = getAlpha5();
				}
				vector<float> Fdash = getFdash(Ftilda, buff_alpha);
				if (!flag) {
					J = lyambda * w + (1 - lyambda) * delta;
					w = getW(Fdash);
					delta = getDelta(Fdash, Ftilda);
					alpha = buff_alpha;
					flag = 1;
				}
				else {
					if (lyambda * w + (1 - lyambda) * delta < J) {
						J = lyambda * w + (1 - lyambda) * delta;
						w = getW(Fdash);
						delta = getDelta(Fdash, Ftilda);
						alpha = buff_alpha;
					}
				}
			}
			print(lyambda, J, alpha, w, delta);
			if (abs(w) + abs(delta) < abs(best_lyambda_w) + abs(best_lyambda_delta)) {
				best_lyambda = lyambda;
				best_lyambda_J = J;
				best_lyambda_w = w;
				best_lyambda_delta = delta;
				best_alpha = alpha;
			}
		}
		cout << endl;
		print_best(best_lyambda, best_lyambda_J, best_lyambda_w, best_lyambda_delta);
		r = 5;
	}
	return 0;
}