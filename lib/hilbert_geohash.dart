library hilbert_geohash;

import 'dart:convert';
import 'dart:math';
import 'package:tuple/tuple.dart';

class HilbertGeoHash {
  final Tuple2 latInterval = const Tuple2<double, double>(-90.0, 90.0);
  final Tuple2 lngInterval = const Tuple2<double, double>(-180.0, 180.0);
  String encode(double lng, double lat, {int precision = 15, int bitsPerChar = 4}) {
    if (lng > lngInterval.item1 &&
        lng < lngInterval.item2 &&
        lat > latInterval.item1 &&
        lat < latInterval.item2 &&
        precision > 0) {
      var bits = precision * bitsPerChar;
      var level = bits >> 1;
      var dim = 1 << level;

      /// INLINE CONVERT COORDS TO INT
      if (dim >= 1) {
        var coords = _coordinatesToInt(lat, lng, dim);


        var code = _intCoordinatesToHash(coords.item1, coords.item2, dim);
        if (code < 0) {
          throw Error();
        }
        var codeString = _encodeInt(code, bitsPerChar);
        return codeString;
      }
    } else {
      throw ArgumentError('Invalid arguments, check if lat/long are within earths bounds or if percision is positive.');
    }
  }

  Tuple2<int, int> _coordinatesToInt(lat, lng, dim) {
    assert(dim >= 1);
    var latY = (lat + latInterval.item2) / 180.0 * dim;
    var lngX = (lng + lngInterval.item2) / 360.0 * dim;

    return Tuple2<int, int>(
        min(dim - 1, int.parse(lngX.floor().toInt().toString())),
        min(dim - 1, int.parse(latY.floor().toInt().toString())));
  }

  int _intCoordinatesToHash(x, y, dim) {
    var d = 0;
    var lvl = dim >> 1;
    while (lvl > 0) {
      int rx = ((x & lvl) > 0) ? 1 : 0;
      int ry = ((y & lvl) > 0) ? 1 : 0;

      d += lvl * lvl * ((3 * rx) ^ ry);

      var xy = _rotateXY(lvl, x, y, rx, ry);
      x = xy.item1;
      y = xy.item2;
      lvl >>= 1;
    }
    return d;
  }

  Tuple2<int, int> _rotateXY(lvl, x, y, rx, ry) {
    if (ry == 0 && rx == 1) {
      var tempy = y;
      y = lvl - 1 - x;
      x = lvl - 1 - tempy;
    } else if (ry == 0) {
      var tempx = x;
      x = y;
      y = tempx;
    }
    return Tuple2(x, y);
  }

  _range(int stop, {int start = 0, int step = 1}) {
    if (step == 0) throw Exception("Step cannot be 0");

    return start < stop == step > 0
        ? List<int>.generate(
            ((start - stop) / step).abs().ceil(), (int i) => start + (i * step))
        : [];
  }

  _encodeInt(code, bitsPerChar) {
    if (code < 0) {
      throw ArgumentError('The hashcode cannot be negative.');
    }
    if (bitsPerChar == 4) {
      var encoded = code.toRadixString(16);
      var encodedSubstring = encoded;
      if (encoded.toString().endsWith('L')) {
        var encodedSubstring = encoded.toString();
        encodedSubstring =
            encodedSubstring.substring(0, encodedSubstring.length);
      } else {
        encodedSubstring = String.fromCharCodes(utf8.encode(encoded));
      }
      return encodedSubstring;
    }
    if (bitsPerChar == 2) {
      var _base4 = '0123';
      var codeLen = ((code.bitLength + 1) / 2).floor();
      var res = ['0'];
      var elements = res.map((e) => e * codeLen).toList();
      for (var i in _range(codeLen - 1, start: -1, step: -1)) {
        elements[i] = _base4[code & 3];
        code >>= 2;
      }
      return elements.toString();
    }
  }

  Tuple2<double, double> decode(code, {int bitsPerChar = 4}) {
    if (code.length == 0) {
      return Tuple2(0.0, 0.0);
    }

    var decodedExactly = decodeExactly(code, bitsPerChar: bitsPerChar);

    return Tuple2(decodedExactly.item1, decodedExactly.item2);
  }

  Map<String, String> neighbours(code, bitsPerChar) {
    var decodedExactly = decodeExactly(code, bitsPerChar: bitsPerChar);
    var precision = code.length;
    var lng = decodedExactly.item1;
    var lat = decodedExactly.item2;
    var lngErr = decodedExactly.item3;
    var latErr = decodedExactly.item4;

    var north = lat + 2 * lngErr;
    var south = lat - 2 * latErr;

    var east = lng + 2 * lngErr;
    if (east > 180) {
      east -= 360;
    }
    var west = lng - 2 * lngErr;
    if (west < -180) {
      west += 360;
    }

    Map neighbursDict = {
      'east': encode(east, lat, precision: precision, bitsPerChar: bitsPerChar),
      'west': encode(west, lat, precision: precision, bitsPerChar: bitsPerChar),
      'north': '',
      'north-east': '',
      'north-west': '',
      'south': '',
      'south-east': '',
      'south-west': '',
    };

    if (north <= 90) {
      neighbursDict.update('north',
          (value) => encode(lng, north, precision: precision, bitsPerChar: bitsPerChar));
      neighbursDict.update('north-east',
          (value) => encode(east, north, precision: precision, bitsPerChar: bitsPerChar));
      neighbursDict.update('north-west',
          (value) => encode(west, north, precision: precision, bitsPerChar: bitsPerChar));
    }

    if (south >= -90) {
      neighbursDict.update('south',
          (value) => encode(lng, south, precision: precision, bitsPerChar: bitsPerChar));
      neighbursDict.update('south-east',
          (value) => encode(east, south, precision: precision, bitsPerChar: bitsPerChar));
      neighbursDict.update('south-west',
          (value) => encode(west, south, precision: precision, bitsPerChar: bitsPerChar));
    }

    return neighbursDict;
  }

  /*
  Expects an object with the following attributes for screenVerticesObject:
      'upLeftLat': double,
      'upLeftLon': double,
      'upRightLat': double,
      'upRightLon': double,
      'downLeftLat': double,
      'downLeftLon': double,
      'downRightLat': double,
      'downRightLon': double,

      level is meant to be in terms of the level of a given geohash cluster.
  */
  List<String> propagateRectangle(Map<String,double> screenVerticesObject, int level, {int precision = 15, int bitsPerChar = 4}) {
    List<String> hashList = [];
    var upLeftHash = encode(screenVerticesObject['upLeftLon'], screenVerticesObject['upLeftLat'], precision: level, bitsPerChar: bitsPerChar);
    var upRightHash = encode(screenVerticesObject['upRightLon'], screenVerticesObject['upRightLat'], precision: level, bitsPerChar: bitsPerChar);
    var downLeftHash = encode(screenVerticesObject['downLeftLon'], screenVerticesObject['downLeftLat'], precision: level, bitsPerChar: bitsPerChar);
    var downRightHash = encode(screenVerticesObject['downRightLon'], screenVerticesObject['downRightLat'], precision: level, bitsPerChar: bitsPerChar);

    var movingHash = upLeftHash;
    var rightLimitHash = upRightHash;
    var endHash = downRightHash;

    // This means the entire screen is in one hash
    if (movingHash == endHash) {
      return hashList;
    }

    var rowHash = movingHash;
    while(movingHash != endHash) { 
      if (movingHash == rightLimitHash) {
        hashList.add(movingHash);
        rowHash = neighbours(rowHash, bitsPerChar)['south'];
        rightLimitHash = neighbours(movingHash,bitsPerChar)['south'];
        movingHash = rowHash;
      }
      while(movingHash != rightLimitHash) {
        var rightNeighbour = neighbours(movingHash,bitsPerChar);
        hashList.add(movingHash);
        movingHash = rightNeighbour['east'];
      }
    } 
    hashList.add(movingHash);
   return hashList;
  }

  Tuple4 decodeExactly(code, {bitsPerChar = 4}) {
    if (bitsPerChar != 2 && bitsPerChar != 4 && bitsPerChar != 6) {
      throw ArgumentError('The base is not supported.');
    }
    if (code.length == 0) {
      return Tuple4(0.0, 0.0, lngInterval.item2, latInterval.item2);
    }
    var bits = code.length * bitsPerChar;
    var level = bits >> 1;
    var dim = 1 << level;

    var codeInt = _decodeInt(code, bitsPerChar);
    var xy = _hashToCoordinates(codeInt, dim);
    var x = xy.item1;
    var y = xy.item2;
    var lnglat = _intToCoordinates(x, y, dim);
    var lngerrlaterr = _lvlError(level);

    return Tuple4(
        lnglat.item1 + lngerrlaterr.item1,
        lnglat.item2 + lngerrlaterr.item2,
        lngerrlaterr.item1,
        lngerrlaterr.item2);
  }

  Tuple2<double, double> _lvlError(lvl) {
    var error = 1 / (1 << lvl);
    return Tuple2(180 * error, 90 * error);
  }

  Tuple2<double, double> _intToCoordinates(x, y, dim) {
    assert(dim >= 1);
    assert(x < dim);
    assert(y < dim);

    var lng = x / dim * 360 - 180;
    var lat = y / dim * 180 - 90;
    return Tuple2(lng, lat);
  }

  Tuple2<int, int> _hashToCoordinates(code, dim) {
    assert(code <= dim * (dim - 1));
    var x = 0;
    var y = 0;
    var lvl = 1;
    while (lvl < dim) {
      var rx = 1 & (code >> 1);
      var ry = 1 & (code ^ rx);
      var xy = _rotateXY(lvl, x, y, rx, ry);
      x = xy.item1;
      y = xy.item2;
      x += lvl * rx;
      y += lvl * ry;
      code >>= 2;
      lvl <<= 1;
    }
    return Tuple2(x, y);
  }

  _decodeInt(tag, bitsPerChar) {
    if (bitsPerChar == 6) {
      return _decodeInt64(tag);
    }
    if (bitsPerChar == 4) {
      return _decodeInt16(tag);
    }
    if (bitsPerChar == 2) {
      return _decodeInt4(tag);
    }
  }

  _decodeInt64(t) {
    if (t.length == 0) {
      return 0;
    }
    return int.parse(t, radix: 64);
  }

  _decodeInt16(t) {
    if (t.length == 0) {
      return 0;
    }
    return int.parse(t, radix: 16);
  }

  _decodeInt4(t) {
    if (t.length == 0) {
      return 0;
    }
    return int.parse(t, radix: 4);
  }
}
