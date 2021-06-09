import { ErrorCorrectionLevels } from "./error_correction_levels.js";
import { QRCodeModes } from "./qr_code_modes.js";
import { BitHandlingUtility } from "./bit_handling_util.js";
import { qrVersionDatabase } from "./qr_version_database.js";
import { characterCountIndicatorLength } from "./character_count_indicator_length.js";
import { ReedSolomon } from "./reed_solomon.js";

const bitHandlingUtilty = new BitHandlingUtility();
export class QRCode {
  constructor(data, error_correction_level) {
    // this data is the string to be encoded
    this.data = data;
    this.error_correction_level = error_correction_level;
    this.ECC_CODEWORDS_PER_BLOCK = [
      // Version: (note that index 0 is for padding, and is set to an illegal value)
      //0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40    Error correction level
      [
        -1, 7, 10, 15, 20, 26, 18, 20, 24, 30, 18, 20, 24, 26, 30, 22, 24, 28,
        30, 28, 28, 28, 28, 30, 30, 26, 28, 30, 30, 30, 30, 30, 30, 30, 30, 30,
        30, 30, 30, 30, 30,
      ], // Low
      [
        -1, 10, 16, 26, 18, 24, 16, 18, 22, 22, 26, 30, 22, 22, 24, 24, 28, 28,
        26, 26, 26, 26, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
        28, 28, 28, 28, 28,
      ], // Medium
      [
        -1, 13, 22, 18, 26, 18, 24, 18, 22, 20, 24, 28, 26, 24, 20, 30, 24, 28,
        28, 26, 30, 28, 30, 30, 30, 30, 28, 30, 30, 30, 30, 30, 30, 30, 30, 30,
        30, 30, 30, 30, 30,
      ], // Quartile
      [
        -1, 17, 28, 22, 16, 22, 28, 26, 26, 24, 28, 24, 28, 22, 24, 24, 30, 28,
        28, 26, 28, 30, 24, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
        30, 30, 30, 30, 30,
      ], // High
    ];
    this.NUM_ERROR_CORRECTION_BLOCKS = [
      // Version: (note that index 0 is for padding, and is set to an illegal value)
      //0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40    Error correction level
      [
        -1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 7, 8, 8, 9, 9,
        10, 12, 12, 12, 13, 14, 15, 16, 17, 18, 19, 19, 20, 21, 22, 24, 25,
      ], // Low
      [
        -1, 1, 1, 1, 2, 2, 4, 4, 4, 5, 5, 5, 8, 9, 9, 10, 10, 11, 13, 14, 16,
        17, 17, 18, 20, 21, 23, 25, 26, 28, 29, 31, 33, 35, 37, 38, 40, 43, 45,
        47, 49,
      ], // Medium
      [
        -1, 1, 1, 2, 2, 4, 4, 6, 6, 8, 8, 8, 10, 12, 16, 12, 17, 16, 18, 21, 20,
        23, 23, 25, 27, 29, 34, 34, 35, 38, 40, 43, 45, 48, 51, 53, 56, 59, 62,
        65, 68,
      ], // Quartile
      [
        -1, 1, 1, 2, 4, 4, 4, 5, 6, 8, 8, 11, 11, 16, 16, 18, 16, 19, 21, 25,
        25, 25, 34, 30, 32, 35, 37, 40, 42, 45, 48, 51, 54, 57, 60, 63, 66, 70,
        74, 77, 81,
      ], // High
    ];
  }

  getQRCodeMode() {
    this.mode = "NUMERIC";
  }

  encodeDataToBitStream(mode, data) {
    //data string will be converted to bit string here
    //numeric mode
    if (mode == "NUMERIC") {
      let resultantCompleteBitString = "";
      let dataLength = data.length;
      for (let i = 0; i < dataLength; i += 3) {
        let currentNumber = parseInt(
          data.substring(i, i + 3 < dataLength ? i + 3 : dataLength)
        );
        resultantCompleteBitString =
          resultantCompleteBitString +
          bitHandlingUtilty.getEncodedBitStringFromNumber(
            currentNumber,
            i + 3 < dataLength
              ? 10
              : currentNumber.toString().length == 1
              ? 4
              : 7,
            true
          );
      }
      this.basic_data_bit_stream = resultantCompleteBitString;
    }
    if (mode == 2) {
      //Alphanumeric code
    }
    if (mode == 3) {
      //Byte Mode
    }
  }

  fitDataToVersion() {
    let version = -1;
    for (version in qrVersionDatabase) {
      let information = qrVersionDatabase[version];
      if (
        parseInt(
          information[this.error_correction_level.toString()]["data_capacity"][
            this.mode
          ]
        ) >= this.data.length
      ) {
        this.version = parseInt(version);
        break;
      }
    }
  }

  getCharacterCountLength() {
    for (let key in characterCountIndicatorLength)
      if (this.version <= parseInt(key))
        return characterCountIndicatorLength[key][this.mode];
  }

  addModeAndCountBits() {
    this.final_data_bit_stream =
      QRCodeModes[this.mode] +
      bitHandlingUtilty.getEncodedBitStringFromNumber(
        this.data.length,
        this.getCharacterCountLength(),
        true
      ) +
      this.basic_data_bit_stream;
  }

  addTerminatorAndBitPadding() {
    let terminatorLength =
      qrVersionDatabase[this.version.toString()][this.error_correction_level][
        "codewords"
      ] *
        8 -
      this.final_data_bit_stream.length;
    if (terminatorLength < 4) {
      for (let i = 0; i < terminatorLength; i++) {
        this.final_data_bit_stream = this.final_data_bit_stream + "0";
      }
    } else {
      this.final_data_bit_stream = this.final_data_bit_stream + "0000";
      let bitPaddingLength = 8 - (this.final_data_bit_stream.length % 8);
      for (let i = 0; i < bitPaddingLength; i++) {
        this.final_data_bit_stream = this.final_data_bit_stream + "0";
      }
    }
  }

  addBytePadding() {
    let information =
      qrVersionDatabase[this.version.toString()][this.error_correction_level];
    let numberOfBytesTillNow = Math.floor(
      this.final_data_bit_stream.length / 8
    );
    let remainingBytes =
      parseInt(information["codewords"]) - numberOfBytesTillNow;
    for (let i = 0; i < remainingBytes / 2; i++) {
      this.final_data_bit_stream =
        this.final_data_bit_stream + "1110110000010001";
    }
    if (remainingBytes % 2 != 0)
      this.final_data_bit_stream = this.final_data_bit_stream + "11101100";
  }

  convertStringBitstreamToArrayBytes() {
    let arr = [];
    for (let i = 0; i < this.final_data_bit_stream.length; i += 8) {
      let currentByte = this.final_data_bit_stream.substring(i, i + 8);
      var digit = parseInt(currentByte, 2);
      arr.push(digit);
    }
    this.finalByteArray = arr;
  }

  addEccAndInterleave() {
    this.finalByteArray = ReedSolomon.addEccAndInterleave(
      this.finalByteArray,
      this
    );
  }
}

let qr = new QRCode("123456789123456789123456789", "Q");
qr.getQRCodeMode();
qr.encodeDataToBitStream(qr.mode, qr.data);
qr.fitDataToVersion();
qr.addModeAndCountBits();
qr.addTerminatorAndBitPadding();
qr.addBytePadding();
qr.convertStringBitstreamToArrayBytes();
qr.addEccAndInterleave();
