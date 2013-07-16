/**
 * Bindable keycodes
 *
 * @property KEY
 * @type Object
 */
gf.input.KEY = {
    BACKSPACE: 8,
    TAB: 9,
    ENTER: 13,
    SHIFT: 16,
    CTRL: 17,
    ALT: 18,
    PAUSE: 19,
    ESC: 27,
    SPACE: 32,
    PAGE_UP: 33,
    PAGE_DOWN: 34,
    END: 35,
    HOME: 36,
    LEFT: 37,
    UP: 38,
    RIGHT: 39,
    DOWN: 40,
    INSERT: 45,
    DELETE: 46,
    NUM0: 48,
    NUM1: 49,
    NUM2: 50,
    NUM3: 51,
    NUM4: 52,
    NUM5: 53,
    NUM6: 54,
    NUM7: 55,
    NUM8: 56,
    NUM9: 57,
    PLUS: 61,
    A : 65,
    B : 66,
    C : 67,
    D : 68,
    E : 69,
    F : 70,
    G : 71,
    H : 72,
    I : 73,
    J : 74,
    K : 75,
    L : 76,
    M : 77,
    N : 78,
    O : 79,
    P : 80,
    Q : 81,
    R : 82,
    S : 83,
    T : 84,
    U : 85,
    V : 86,
    W : 87,
    X : 88,
    Y : 89,
    Z : 90,
    NUMPAD0: 96,
    NUMPAD1: 97,
    NUMPAD2: 98,
    NUMPAD3: 99,
    NUMPAD4: 100,
    NUMPAD5: 101,
    NUMPAD6: 102,
    NUMPAD7: 103,
    NUMPAD8: 104,
    NUMPAD9: 105,
    NUMPAD_STAR: 106,
    NUMPAD_PLUS: 107,
    NUMPAD_MINUS: 109,
    NUMPAD_DOT: 110,
    NUMPAD_SLASH: 111,
    F1: 112,
    F2: 113,
    F3: 114,
    F4: 115,
    MINUS: 173,
    TILDE: 192
};

/**
 * Controls keyboard input
 *
 * @class Keyboard
 * @namespace input
 * @constructor
 * @param view {DOMElement} The DOMElement to bind input events to
 */
gf.input.Keyboard = function(view) {
    gf.input.Input.call(this, view);

    /**
     * The current sequence of keys that have been pressed
     *
     * @property sequence
     * @type Array<Number>
     * @readOnly
     */
    this.sequence = [];

    /**
     * The amount of time it takes for the sequence to clear out, in ms
     *
     * @property sequenceTimeout
     * @type Number
     * @default 500
     */
    this.sequenceTimeout = 500;

    /**
     * The timeout ID for the wait to clear the input sequence
     *
     * @property _clearSq
     * @type Number
     * @private 
     */
    this._clearSq = null;

    document.addEventListener('keydown', this.onKeyDown.bind(this), false);
    document.addEventListener('keyup', this.onKeyUp.bind(this), false);
};

gf.inherits(gf.input.Keyboard, gf.input.Input, {
    //on keydown event set gf.controls keycode's action as active
    //and call any registered callbacks
    onKeyDown: function(e, override) {
        if(e.target === this.view)
            return this.modifyKey(e, override || e.keyCode || e.which, true);
    },
    onKeyUp: function(e, override) {
        if(e.target === this.view)
            return this.modifyKey(e, override || e.keyCode || e.which, false);
    },
    modifyKey: function(e, key, down) {
        //emit single key event
        this.emit(key, {
            originalEvent: e,
            down: down
        });

        //when pressed is when we process a key for a sequence
        if(down) {
            //update the key sequence
            this.sequence.push(key);

            //process current sequence
            this.emit(this.sequence.toString(), {
                originalEvent: e,
                down: down
            });

            //set timeout to clear sequence
            clearTimeout(this._clearSq);
            this._clearSq = setTimeout(this._clearSequence.bind(this), this.sequenceTimeout);
        }
    },
    _clearSequence: function() {
        this.sequence.length = 0;
    }
});