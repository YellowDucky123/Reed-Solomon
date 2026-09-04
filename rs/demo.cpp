#include "reed_solomon.hpp"
#include <iostream>
#include <string>

// Small end-to-end demonstration: encode a message, corrupt a few symbols of
// the codeword, then decode and recover the original.
int main() {
	std::string message = "REEDSOLOMON";
	int parity = 6; // corrects up to parity/2 = 3 symbol errors

	ReedSolomon rs((int)message.size(), parity);
	std::vector<Fp> code = rs.encode(message);

	std::cout << "message : " << message << "\n";
	std::cout << "k=" << rs.k << " parity=" << rs.parity << " n=" << rs.n
	          << " correctable=" << rs.correctable() << "\n";

	std::cout << "codeword:";
	for (const Fp& c : code) std::cout << " " << c.v;
	std::cout << "\n";

	// Corrupt 3 symbols.
	code[1] = code[1] + Fp(123);
	code[4] = code[4] + Fp(9999);
	code[10] = code[10] + Fp(4);
	std::cout << "corrupted 3 symbols (positions 1, 4, 10)\n";

	ReedSolomon::DecodeString d = rs.decode(code);
	std::cout << "decoded : " << (d.ok ? d.message : "<decode failed>")
	          << (d.ok && d.message == message ? "  [recovered]" : "") << "\n";
	return d.ok && d.message == message ? 0 : 1;
}
