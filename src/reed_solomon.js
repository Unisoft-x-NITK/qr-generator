"use strict";
import { ErrorCorrectionLevels } from "./error_correction_levels.js";

    export var ReedSolomon = /** @class */ (function () {
        function ReedSolomon() {
        }
        ReedSolomon.getNumRawDataModules = function (ver) {
            if (ver < 1 || ver > 40)
                throw "Version number out of range";
            var result = (16 * ver + 128) * ver + 64;
            if (ver >= 2) {
                var numAlign = Math.floor(ver / 7) + 2;
                result -= (25 * numAlign - 10) * numAlign - 55;
                if (ver >= 7)
                    result -= 36;
            }
            if (!(208 <= result && result <= 29648))
                throw "Assertion error";
            return result;
        };
        // Returns the number of 8-bit data (i.e. not error correction) codewords contained in any
        // QR Code of the given version number and error correction level, with remainder bits discarded.
        // This stateless pure function could be implemented as a (40*4)-cell lookup table.
        ReedSolomon.getNumDataCodewords = function (ver, QrCode) {
            var ecl = QrCode.error_correction_level;
            return Math.floor(ReedSolomon.getNumRawDataModules(ver) / 8) -
                QrCode.ECC_CODEWORDS_PER_BLOCK[ErrorCorrectionLevels[ecl]][ver] *
                    QrCode.NUM_ERROR_CORRECTION_BLOCKS[ErrorCorrectionLevels[ecl]][ver];
        };
        // Returns a Reed-Solomon ECC generator polynomial for the given degree. This could be
        // implemented as a lookup table over all possible parameter values, instead of as an algorithm.
        ReedSolomon.reedSolomonComputeDivisor = function (degree) {
            if (degree < 1 || degree > 255)
                throw "Degree out of range";
            // Polynomial coefficients are stored from highest to lowest power, excluding the leading term which is always 1.
            // For example the polynomial x^3 + 255x^2 + 8x + 93 is stored as the uint8 array [255, 8, 93].
            var result = [];
            for (var i = 0; i < degree - 1; i++)
                result.push(0);
            result.push(1); // Start off with the monomial x^0
            // Compute the product polynomial (x - r^0) * (x - r^1) * (x - r^2) * ... * (x - r^{degree-1}),
            // and drop the highest monomial term which is always 1x^degree.
            // Note that r = 0x02, which is a generator element of this field GF(2^8/0x11D).
            var root = 1;
            for (var i = 0; i < degree; i++) {
                // Multiply the current product by (x - r^i)
                for (var j = 0; j < result.length; j++) {
                    result[j] = ReedSolomon.reedSolomonMultiply(result[j], root);
                    if (j + 1 < result.length)
                        result[j] ^= result[j + 1];
                }
                root = ReedSolomon.reedSolomonMultiply(root, 0x02);
            }
            return result;
        };
        // Returns the Reed-Solomon error correction codeword for the given data and divisor polynomials.
        ReedSolomon.reedSolomonComputeRemainder = function (data, divisor) {
            var result = divisor.map(function (_) { return 0; });
            var _loop_1 = function (b) {
                var factor = b ^ result.shift();
                result.push(0);
                divisor.forEach(function (coef, i) {
                    return result[i] ^= ReedSolomon.reedSolomonMultiply(coef, factor);
                });
            };
            for (var _i = 0, data_1 = data; _i < data_1.length; _i++) {
                var b = data_1[_i];
                _loop_1(b);
            }
            return result;
        };
        // Returns the product of the two given field elements modulo GF(2^8/0x11D). The arguments and result
        // are unsigned 8-bit integers. This could be implemented as a lookup table of 256*256 entries of uint8.
        ReedSolomon.reedSolomonMultiply = function (x, y) {
            if (x >>> 8 != 0 || y >>> 8 != 0)
                throw "Byte out of range";
            // Russian peasant multiplication
            var z = 0;
            for (var i = 7; i >= 0; i--) {
                z = (z << 1) ^ ((z >>> 7) * 0x11D);
                z ^= ((y >>> i) & 1) * x;
            }
            if (z >>> 8 != 0)
                throw "Assertion error";
            return z;
        };
        ReedSolomon.addEccAndInterleave = function (data, QrCode) {
            var ver = QrCode.version;
            var ecl = QrCode.error_correction_level;
            if (data.length != ReedSolomon.getNumDataCodewords(ver, QrCode))
                throw "Invalid argument";
            // Calculate parameter numbers
            var numBlocks = QrCode.NUM_ERROR_CORRECTION_BLOCKS[ErrorCorrectionLevels[ecl]][ver];
            var blockEccLen = QrCode.ECC_CODEWORDS_PER_BLOCK[ErrorCorrectionLevels[ecl]][ver];
            var rawCodewords = Math.floor(ReedSolomon.getNumRawDataModules(ver) / 8);
            var numShortBlocks = numBlocks - rawCodewords % numBlocks;
            var shortBlockLen = Math.floor(rawCodewords / numBlocks);
            // Split data into blocks and append ECC to each block
            var blocks = [];
            var rsDiv = ReedSolomon.reedSolomonComputeDivisor(blockEccLen);
            for (var i = 0, k = 0; i < numBlocks; i++) {
                var dat = data.slice(k, k + shortBlockLen - blockEccLen + (i < numShortBlocks ? 0 : 1));
                k += dat.length;
                var ecc = ReedSolomon.reedSolomonComputeRemainder(dat, rsDiv);
                if (i < numShortBlocks)
                    dat.push(0);
                blocks.push(dat.concat(ecc));
            }
            // Interleave (not concatenate) the bytes from every block into a single sequence
            var result = [];
            var _loop_2 = function (i) {
                blocks.forEach(function (block, j) {
                    // Skip the padding byte in short blocks
                    if (i != shortBlockLen - blockEccLen || j >= numShortBlocks)
                        result.push(block[i]);
                });
            };
            for (var i = 0; i < blocks[0].length; i++) {
                _loop_2(i);
            }
            if (result.length != rawCodewords)
                throw "Assertion error";
            return result;
        };
        return ReedSolomon;
    }());