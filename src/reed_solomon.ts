import { ErrorCorrectionLevels } from './error_correction_levels.js';
import {QRCode} from './qr.js';
namespace reedsolomonnamespace {
    type bit = number;
    type byte = number;
    type int = number;
    export class ReedSolomon {
        private static getNumRawDataModules(ver): int {
            if (ver < 1 || ver > 40)
                throw "Version number out of range";
            let result: int = (16 * ver + 128) * ver + 64;
            if (ver >= 2) {
                const numAlign: int = Math.floor(ver / 7) + 2;
                result -= (25 * numAlign - 10) * numAlign - 55;
                if (ver >= 7)
                    result -= 36;
            }
            if (!(208 <= result && result <= 29648))
                throw "Assertion error";
            return result;
        }


        // Returns the number of 8-bit data (i.e. not error correction) codewords contained in any
        // QR Code of the given version number and error correction level, with remainder bits discarded.
        // This stateless pure function could be implemented as a (40*4)-cell lookup table.
        private static getNumDataCodewords(ver: int, QrCode:QRCode): int {
            let ecl = QrCode.error_correction_level;
            return Math.floor(ReedSolomon.getNumRawDataModules(ver) / 8) -
                QrCode.ECC_CODEWORDS_PER_BLOCK[ErrorCorrectionLevels[ecl]][ver] *
                QrCode.NUM_ERROR_CORRECTION_BLOCKS[ErrorCorrectionLevels[ecl]][ver];
        }


        // Returns a Reed-Solomon ECC generator polynomial for the given degree. This could be
        // implemented as a lookup table over all possible parameter values, instead of as an algorithm.
        private static reedSolomonComputeDivisor(degree: int): Array<byte> {
            if (degree < 1 || degree > 255)
                throw "Degree out of range";
            // Polynomial coefficients are stored from highest to lowest power, excluding the leading term which is always 1.
            // For example the polynomial x^3 + 255x^2 + 8x + 93 is stored as the uint8 array [255, 8, 93].
            let result: Array<byte> = [];
            for (let i = 0; i < degree - 1; i++)
                result.push(0);
            result.push(1);  // Start off with the monomial x^0

            // Compute the product polynomial (x - r^0) * (x - r^1) * (x - r^2) * ... * (x - r^{degree-1}),
            // and drop the highest monomial term which is always 1x^degree.
            // Note that r = 0x02, which is a generator element of this field GF(2^8/0x11D).
            let root = 1;
            for (let i = 0; i < degree; i++) {
                // Multiply the current product by (x - r^i)
                for (let j = 0; j < result.length; j++) {
                    result[j] = ReedSolomon.reedSolomonMultiply(result[j], root);
                    if (j + 1 < result.length)
                        result[j] ^= result[j + 1];
                }
                root = ReedSolomon.reedSolomonMultiply(root, 0x02);
            }
            return result;
        }


        // Returns the Reed-Solomon error correction codeword for the given data and divisor polynomials.
        private static reedSolomonComputeRemainder(data: Array<byte>, divisor: Array<byte>): Array<byte> {
            let result: Array<byte> = divisor.map(_ => 0);
            for (const b of data) {  // Polynomial division
                const factor: byte = b ^ (result.shift() as byte);
                result.push(0);
                divisor.forEach((coef, i) =>
                    result[i] ^= ReedSolomon.reedSolomonMultiply(coef, factor));
            }
            return result;
        }


        // Returns the product of the two given field elements modulo GF(2^8/0x11D). The arguments and result
        // are unsigned 8-bit integers. This could be implemented as a lookup table of 256*256 entries of uint8.
        private static reedSolomonMultiply(x: byte, y: byte): byte {
            if (x >>> 8 != 0 || y >>> 8 != 0)
                throw "Byte out of range";
            // Russian peasant multiplication
            let z: int = 0;
            for (let i = 7; i >= 0; i--) {
                z = (z << 1) ^ ((z >>> 7) * 0x11D);
                z ^= ((y >>> i) & 1) * x;
            }
            if (z >>> 8 != 0)
                throw "Assertion error";
            return z as byte;
        }
        private addEccAndInterleave(data: Array<byte>, QrCode: QRCode): Array<byte> {
            const ver: int = QrCode.version;
            const ecl: string = QrCode.error_correction_level;
            if (data.length != ReedSolomon.getNumDataCodewords(ver, QrCode))
                throw "Invalid argument";

            // Calculate parameter numbers
            const numBlocks: int = QrCode.NUM_ERROR_CORRECTION_BLOCKS[ErrorCorrectionLevels[ecl]][ver];
            const blockEccLen: int = QrCode.ECC_CODEWORDS_PER_BLOCK[ErrorCorrectionLevels[ecl]][ver];
            const rawCodewords: int = Math.floor(ReedSolomon.getNumRawDataModules(ver) / 8);
            const numShortBlocks: int = numBlocks - rawCodewords % numBlocks;
            const shortBlockLen: int = Math.floor(rawCodewords / numBlocks);

            // Split data into blocks and append ECC to each block
            let blocks: Array<Array<byte>> = [];
            const rsDiv: Array<byte> = ReedSolomon.reedSolomonComputeDivisor(blockEccLen);
            for (let i = 0, k = 0; i < numBlocks; i++) {
                let dat: Array<byte> = data.slice(k, k + shortBlockLen - blockEccLen + (i < numShortBlocks ? 0 : 1));
                k += dat.length;
                const ecc: Array<byte> = ReedSolomon.reedSolomonComputeRemainder(dat, rsDiv);
                if (i < numShortBlocks)
                    dat.push(0);
                blocks.push(dat.concat(ecc));
            }

            // Interleave (not concatenate) the bytes from every block into a single sequence
            let result: Array<byte> = [];
            for (let i = 0; i < blocks[0].length; i++) {
                blocks.forEach((block, j) => {
                    // Skip the padding byte in short blocks
                    if (i != shortBlockLen - blockEccLen || j >= numShortBlocks)
                        result.push(block[i]);
                });
            }
            if (result.length != rawCodewords)
                throw "Assertion error";
            return result;
        }
    }
}
